// lib/wqrs.dart
import 'dart:math';
import 'dart:io';

class WQRS {
  // Константы из оригинального C-кода
  static const int BUFLN = 16384;
  static const double EYE_CLS = 0.25;
  static const double MaxQRSw = 0.13;
  static const double NDP = 2.5;
  static const int PWFreqDEF = 60;
  static const int TmDEF = 100;

  // Параметры класса
  double _sampleRate = 250.0;
  int _sig = -1;
  int _Tm = TmDEF;
  int _PWFreq = PWFreqDEF;
  
  // Параметры, вычисляемые на основе частоты дискретизации
  late int _LPn;
  late int _LP2n;
  late int _LTwindow;
  late int _EyeClosing;
  late int _ExpectPeriod;
  late double _lfsc;

  // Состояние ltsamp
  List<int>? _lbuf;
  List<int>? _ebuf;
  int _tt = -1;
  int _Yn = 0;
  int _Yn1 = 0;
  int _Yn2 = 0;
  int _aet = 0;
  int _et = 0;

  // Входной сигнал
  List<double>? _signal;
  double? _gain;

  WQRS({
    double sampleRate = 250.0,
    int signal = -1,
    int threshold = TmDEF,
    int powerFreq = PWFreqDEF,
    double? gain,
  }) {
    _sampleRate = sampleRate;
    _sig = signal;
    _Tm = threshold;
    _PWFreq = powerFreq;
    _gain = gain;
    _initializeParameters();
    reset();
  }

  void _initializeParameters() {
    // Фильтры для подавления сетевой наводки (50 или 60 Гц)
    _LPn = (_sampleRate / _PWFreq).round();
    if (_LPn > 8) _LPn = 8;
    _LP2n = 2 * _LPn;
    
    // Параметры детекции
    _EyeClosing = (_sampleRate * EYE_CLS).round();
    _ExpectPeriod = (_sampleRate * NDP).round();
    _LTwindow = (_sampleRate * MaxQRSw).round();
    
    _lfsc = 1.25 * _gain! * _gain! / _sampleRate;
    
    print('WQRS initialized:');
    print('  Sample rate: $_sampleRate Hz');
    print('  Power line frequency: $_PWFreq Hz');
    print('  LPn: $_LPn, LP2n: $_LP2n');
    print('  lfsc: $_lfsc');
    print('  Eye closing: $_EyeClosing samples');
    print('  LT window: $_LTwindow samples');
  }

  int _sample(int sig, int t) {
    if (_signal == null || t < 0 || t >= _signal!.length) {
      return -1; // WFDB_INVALID_SAMPLE
    }
    // (в оригинале WFDB_Sample это целое число)
    return _signal![t].round();
  }

  bool _sampleValid() {
    return _signal != null && _signal!.isNotEmpty;
  }

  int ltsamp(int t) {
    int dy;

    // Инициализация буферов
    if (_lbuf == null) {
      _lbuf = List<int>.filled(BUFLN, 0);
      _ebuf = List<int>.filled(BUFLN, 0);

      if (_lbuf != null && _ebuf != null) {
        _ebuf![0] = sqrt(_lfsc).round();
        _tt = 1;
        while (_tt < BUFLN) {
          _ebuf![_tt] = _ebuf![0];
          _tt++;
        }

        if (t > BUFLN) {
          _tt = t - BUFLN;
        } else {
          _tt = -1;
        }

        _Yn = 0;
        _Yn1 = 0;
        _Yn2 = 0;
      } else {
        stderr.writeln('WQRS: insufficient memory');
        exit(2);
      }
    }

    if (t < _tt - BUFLN) {
      stderr.writeln('WQRS: ltsamp buffer too short');
      exit(2);
    }

    while (t > _tt) {
      _Yn2 = _Yn1;
      _Yn1 = _Yn;

      int v0 = _sample(_sig, _tt);
      int v1 = _sample(_sig, _tt - _LPn);
      int v2 = _sample(_sig, _tt - _LP2n);

      if (v0 != -1 && v1 != -1 && v2 != -1) {
        _Yn = 2 * _Yn1 - _Yn2 + v0 - 2 * v1 + v2;
      }

      dy = (_Yn - _Yn1) ~/ _LP2n;
      _tt++;

      _et = sqrt(_lfsc + dy * dy).round();
      _ebuf![_tt & (BUFLN - 1)] = _et;

      _aet += _et - _ebuf![(_tt - _LTwindow) & (BUFLN - 1)];
      _lbuf![_tt & (BUFLN - 1)] = _aet;
    }

    return _lbuf![t & (BUFLN - 1)];
  }

  void reset() {
    _lbuf = null;
    _ebuf = null;
    _tt = -1;
    _Yn = 0;
    _Yn1 = 0;
    _Yn2 = 0;
    _aet = 0;
    _et = 0;
    _signal = null;
  }

  // Основной метод обработки
  List<Map<String, dynamic>> process(List<double> signal, {int? startTime, int? endTime}) {
    if (signal.isEmpty) return [];

    reset();
    _signal = signal;

    int from = startTime ?? 0;
    int to = endTime ?? signal.length;

    List<Map<String, dynamic>> detections = [];

    // Вычисляем начальные пороги (усредняем первые 8 секунд)
    int t1 = (8 * _sampleRate).round();
    if (t1 > BUFLN * 0.9) t1 = (BUFLN / 2).round();
    t1 += from;

    double T0 = 0;
    int sampleCount = 0;
    for (int t = from; t < t1 && t < signal.length; t++) {
      T0 += ltsamp(t);
      sampleCount++;
    }
    if (sampleCount > 0) T0 /= sampleCount;
    double Ta = 3 * T0;

    bool learning = true;
    int timer = 0;
    int minutes = 0;
    int next_minute = from + (60 * _sampleRate).round();
    double T1 = 0;

    for (int t = from; t < to && t < signal.length; t++) {
      if (learning) {
        if (t > t1) {
          learning = false;
          T1 = T0;
          t = from;
        } else {
          T1 = 2 * T0;
        }
      }

      if (ltsamp(t) > T1) {
        timer = 0;
        int maxVal = ltsamp(t);
        int minVal = ltsamp(t);
        
        // Поиск максимума и минимума в окне
        for (int tt = t + 1; tt < t + (_EyeClosing / 2).round() && tt < signal.length; tt++) {
          int val = ltsamp(tt);
          if (val > maxVal) maxVal = val;
        }
        for (int tt = t - 1; tt > t - (_EyeClosing / 2).round() && tt >= 0; tt--) {
          int val = ltsamp(tt);
          if (val < minVal) minVal = val;
        }
        
        if (maxVal > minVal + 10) {
          // Находим начало QRS (PQ junction)
          int onset = (maxVal / 100).round() + 2;
          int tpq = t - 5;
          for (int tt = t; tt > t - (_EyeClosing / 2).round() && tt >= 4; tt--) {
            if (ltsamp(tt) - ltsamp(tt - 1) < onset &&
                ltsamp(tt - 1) - ltsamp(tt - 2) < onset &&
                ltsamp(tt - 2) - ltsamp(tt - 3) < onset &&
                ltsamp(tt - 3) - ltsamp(tt - 4) < onset) {
              tpq = tt - _LP2n;
              break;
            }
          }

          if (!learning) {
            _sample(_sig, tpq);
            if (!_sampleValid()) break;
            
            detections.add({
              'time': tpq,
              'type': 'QRS_ONSET',
              'value': ltsamp(t)
            });
          }

          // Обновляем пороги
          Ta += (maxVal - Ta) / 10;
          T1 = Ta / 3;

          // Блокировка на время eye-closing
          t += _EyeClosing;
        }
      } else if (!learning) {
        timer++;
        if (timer > _ExpectPeriod && Ta > _Tm) {
          Ta--;
          T1 = Ta / 3;
        }
      }

      // Отслеживание прогресса
      if (t >= next_minute) {
        next_minute += (60 * _sampleRate).round();
        minutes++;
        if (minutes % 10 == 0) {
          print('Обработано $minutes минут');
        }
      }
    }

    return detections;
  }

  // Вспомогательный метод для получения значений ltsamp
  List<int> getLTSampValues(List<double> signal) {
    _signal = signal;
    List<int> result = [];
    for (int i = 0; i < signal.length; i++) {
      result.add(ltsamp(i));
    }
    return result;
  }

  void debugLTSamp(List<double> signal, int start, int length) {
    _signal = signal;
    print('Диагностика ltsamp для первых $length сэмплов:');
    print('sample\tltsamp');
    for (int i = start; i < start + length && i < signal.length; i++) {
      int val = ltsamp(i);
      if (i < start + 10 || i % 100 == 0) {
        print('$i\t$val');
      }
    }
  }

}