// lib/wqrs_detector.dart
import 'dart:math';
import 'dart:io';
/// Класс для обнаружения QRS комплексов на основе алгоритма wqrs
/// Работает с сырыми АЦП-данными (оригинальный режим)
class WQRSDetector {
  // Параметры алгоритма
  final double _sampleRate;
  final double _lfsc;
  final int _lpn;
  final int _lp2n;
  final int _eyeClosingSamples;
  final int _expectPeriod;
  final int _minThreshold;
  final int _ltWindow;
  
  // Буферы
  static const int _BUFFER_SIZE = 16384;
  final List<int> _buffer = List.filled(_BUFFER_SIZE, 0);
  final List<int> _ebuf = List.filled(_BUFFER_SIZE, 0);
  final List<int> _lbuf = List.filled(_BUFFER_SIZE, 0);
  
  // Состояния фильтра
  int _yn = 0;
  int _yn1 = 0;
  int _yn2 = 0;
  int _tt = -1;
  int _aet = 0;
  
  // Пороги
  double _t0 = 0;  // Изменено на double для точности
  double _ta = 0;  // Изменено на double для точности
  int _t1 = 0;
  
  // Состояние обучения
  bool _isLearning = true;
  int _t1_learning = 0;
  
  // Состояние детекции
  int _lastDetection = -1000;
  int _timer = 0;
  bool _inRefractory = false;
  
  // Результаты
  final List<int> _detections = [];
  final List<int> _jPoints = [];
  
  // Параметры
  final double _gain;
  final bool _debugMode;
  final bool _inputIsRaw;
  int _samplesProcessed = 0;
  
  // Для отслеживания прогресса
  int _minutes = 0;
  int _nextMinute = 0;
  
  WQRSDetector({
    required double sampleRate,
    double gain = 200.0,
    int powerLineFreq = 50,
    int minThreshold = 100,
    double eyeClosingPeriod = 0.25,
    double maxQRSWidth = 0.13,
    double ndp = 2.5,
    bool debugMode = false,
    bool inputIsRaw = true,
  }) : _sampleRate = sampleRate,
       _gain = gain,
       _minThreshold = minThreshold,
       _lpn = (sampleRate / powerLineFreq).floor().clamp(1, 8),
       _lp2n = (2 * (sampleRate / powerLineFreq).floor()).clamp(2, 16),
       _eyeClosingSamples = (eyeClosingPeriod * sampleRate).round(),
       _expectPeriod = (ndp * sampleRate).round(),
       _ltWindow = (maxQRSWidth * sampleRate).round().clamp(1, 256),
       _lfsc = 1.25 * gain * gain / sampleRate,
       _debugMode = debugMode,
       _inputIsRaw = inputIsRaw {
    
    // Инициализация буфера ebuf как в оригинале
    int initialValue = sqrt(_lfsc).round();
    for (int i = 0; i < _BUFFER_SIZE; i++) {
      _ebuf[i] = initialValue;
    }
    
    if (_debugMode) {
      print('=== WQRS Detector Initialized ===');
      print('Sample Rate: $_sampleRate Hz');
      print('Gain: $_gain ADC/ед.');
      print('Input mode: ${_inputIsRaw ? "RAW ADC" : "mV"}');
      print('LPn: $_lpn, LP2n: $_lp2n');
      print('Eye Closing Samples: $_eyeClosingSamples');
      print('LT Window: $_ltWindow');
      print('LFSC: $_lfsc');
      print('Min Threshold: $_minThreshold');
      print('================================');
    }
  }
  
  void reset() {
    _buffer.fillRange(0, _BUFFER_SIZE, 0);
    int initialValue = sqrt(_lfsc).round();
    _ebuf.fillRange(0, _BUFFER_SIZE, initialValue);
    _lbuf.fillRange(0, _BUFFER_SIZE, 0);
    
    _yn = 0;
    _yn1 = 0;
    _yn2 = 0;
    _tt = -1;
    _aet = 0;
    
    _t0 = 0;
    _ta = 0;
    _t1 = 0;
    _t1_learning = 0;
    
    _isLearning = true;
    _lastDetection = -1000;
    _timer = 0;
    _inRefractory = false;
    
    _detections.clear();
    _jPoints.clear();
    _samplesProcessed = 0;
    _minutes = 0;
    _nextMinute = (8 * _sampleRate).round(); // 8 секунд как в оригинале
  }
  
  int _toAdc(dynamic value) {
    if (_inputIsRaw) {
      return value is int ? value : value.round();
    } else {
      return (value * _gain).round();
    }
  }
  
  int _ltsamp(int t, List<int> signal) {
    int dy;
    int v0, v1, v2;
    int et;
    
    if (_tt == -1) {
      int initialValue = sqrt(_lfsc).round();
      for (int i = 0; i < _BUFFER_SIZE; i++) {
        _ebuf[i] = initialValue;
      }
      _tt = (t > _BUFFER_SIZE) ? t - _BUFFER_SIZE : -1;
      _yn = _yn1 = _yn2 = 0;
      _aet = 0;
    }
    
    if (t < _tt - _BUFFER_SIZE) {
      if (_debugMode) print('ERROR: ltsamp buffer too short at t=$t, tt=$_tt');
      return 0;
    }
    
    while (t > _tt) {
      _yn2 = _yn1;
      _yn1 = _yn;
      
      v0 = (_tt >= 0 && _tt < signal.length) ? signal[_tt] : 0;
      v1 = (_tt - _lpn >= 0 && _tt - _lpn < signal.length) ? signal[_tt - _lpn] : 0;
      v2 = (_tt - _lp2n >= 0 && _tt - _lp2n < signal.length) ? signal[_tt - _lp2n] : 0;
      
      // Фильтр низких частот (оригинальный)
      _yn = 2 * _yn1 - _yn2 + v0 - 2 * v1 + v2;
      
      // Производная
      dy = (_yn - _yn1) ~/ _lp2n;
      
      // Преобразование длины
      _tt++;
      double sqrtVal = sqrt(_lfsc + dy * dy);
      et = sqrtVal.round();
      
      _ebuf[_tt & (_BUFFER_SIZE - 1)] = et;
      _aet += et - _ebuf[(_tt - _ltWindow) & (_BUFFER_SIZE - 1)];
      _lbuf[_tt & (_BUFFER_SIZE - 1)] = _aet;
    }
    
    return _lbuf[t & (_BUFFER_SIZE - 1)];
  }
  
  // Новый метод для обработки одного семпла, возвращает true если обнаружен QRS
  bool _processSample(int sample, List<int> signal, int t, bool isLearningPhase) {
    int lt = _ltsamp(t, signal);
    
    // === ФАЗА ОБУЧЕНИЯ ===
    if (isLearningPhase) {
      // В оригинале T1 = 2*T0 во время обучения, но мы не детектим
      return false;
    }
    
    // === ФАЗА ДЕТЕКЦИИ ===
    // Проверка рефрактерного периода
    if (_inRefractory) {
      // В оригинале просто пропускаем семплы через t += EyeClosing
      // Здесь мы эмулируем это через проверку
      return false;
    }
    
    // Отладочный вывод
    if (_debugMode && t % 5000 == 0) {
      print('Sample $t: LT=$lt, T1=$_t1, detections=${_detections.length}');
    }
    
    // Сравнение с порогом
    if (lt > _t1) {
      _timer = 0;
      
      // Поиск максимума в окне вперед
      int maxVal = lt;
      int halfWindow = _eyeClosingSamples ~/ 2;
      
      for (int tt = t + 1; tt < t + halfWindow && tt < signal.length; tt++) {
        int val = _ltsamp(tt, signal);
        if (val > maxVal) maxVal = val;
      }
      
      // Поиск минимума в окне назад
      int minVal = lt;
      for (int tt = t - 1; tt > t - halfWindow && tt >= 0; tt--) {
        int val = _ltsamp(tt, signal);
        if (val < minVal) minVal = val;
      }
      
      // Проверка наличия QRS
      if (maxVal > minVal + 10) {
        // Поиск начала QRS (PQ junction)
        int onset = (maxVal ~/ 100) + 2;
        int tpq = t - 5;
        
        for (int tt = t; tt > t - halfWindow && tt >= 4; tt--) {
          int val0 = _ltsamp(tt, signal);
          int val1 = _ltsamp(tt - 1, signal);
          int val2 = _ltsamp(tt - 2, signal);
          int val3 = _ltsamp(tt - 3, signal);
          int val4 = _ltsamp(tt - 4, signal);
          
          if (val0 - val1 < onset &&
              val1 - val2 < onset &&
              val2 - val3 < onset &&
              val3 - val4 < onset) {
            tpq = tt - _lp2n;
            break;
          }
        }
        
        // Регистрация обнаружения (только если не в режиме обучения)
        if (tpq >= 0 && tpq < signal.length) {
          if (_detections.isEmpty || tpq - _detections.last > _eyeClosingSamples ~/ 2) {
            _detections.add(tpq);
            _lastDetection = t;
            _inRefractory = true;
            
            if (_debugMode && _detections.length % 100 == 0) {
              print('>>> QRS DETECTED #${_detections.length} at index $tpq');
            }
            
            // Поиск J-точки
            int jPoint = t + 5;
            for (int tt = t; tt < t + halfWindow && tt < signal.length; tt++) {
              int val = _ltsamp(tt, signal);
              if (val > maxVal - (maxVal ~/ 10)) {
                jPoint = tt;
                break;
              }
            }
            if (jPoint < signal.length) {
              _jPoints.add(jPoint);
            }
            
            // Адаптация порогов
            _ta += (maxVal - _ta) / 10;
            _t1 = (_ta / 3).round();
          }
        }
        
        // В оригинале здесь t += EyeClosing (пропуск семплов)
        // Мы эмулируем это через флаг рефрактерности и возвращаем true
        return true;
      }
    } else {
      // Уменьшение порога
      _timer++;
      if (_timer > _expectPeriod && _ta > _minThreshold) {
        _ta--;
        _t1 = (_ta / 3).round();
        _timer = 0;
        
        if (_debugMode && _ta % 100 == 0) {
          print('Threshold lowered: Ta=$_ta, T1=$_t1');
        }
      }
    }
    
    return false;
  }
  
  List<int> detect(List<int> signal) {
    if (_debugMode) {
      print('=== Starting Detection ===');
      print('Signal length: ${signal.length} samples');
      print('Input mode: ${_inputIsRaw ? "RAW ADC" : "mV"}');
      
      if (signal.isNotEmpty) {
        int minSignal = signal.reduce((a, b) => a < b ? a : b);
        int maxSignal = signal.reduce((a, b) => a > b ? a : b);
        print('Signal range: $minSignal to $maxSignal ADC');
        print('Signal amplitude: ${maxSignal - minSignal} ADC');
        
        if (_inputIsRaw) {
          if (maxSignal - minSignal < 100) {
            print('WARNING: Signal amplitude is very small!');
            print('  Expected ADC amplitude for ECG: ~2000-5000');
          }
        }
      }
    }
    
    reset();
    
    // === ОБУЧЕНИЕ (как в оригинале) ===
    if (_debugMode) print('Phase 1: Learning...');
    
    // Вычисляем t1 = 8 секунд, но не более половины буфера
    int t1 = (8.0 * _sampleRate).round();
    if (t1 > _BUFFER_SIZE * 0.9) {
      t1 = _BUFFER_SIZE ~/ 2;
    }
    t1 = t1.clamp(1, signal.length);
    
    // Накопление T0 за первые 8 секунд
    double T0 = 0;
    int samplesCount = 0;
    
    for (int t = 0; t < t1 && t < signal.length; t++) {
      int lt = _ltsamp(t, signal);
      T0 += lt;
      samplesCount++;
    }
    
    if (samplesCount == 0) {
      if (_debugMode) print('ERROR: Signal too short for learning');
      return [];
    }
    
    T0 /= samplesCount;  // T0 = среднее значение
    
    if (_debugMode) {
      print('Learning complete: T0=$T0, samples=$samplesCount');
    }
    
    // === ПЕРЕЗАПУСК (как в оригинале) ===
    if (_debugMode) print('Phase 2: Restarting and detecting...');
    
    // Полный сброс состояния
    _yn = 0;
    _yn1 = 0;
    _yn2 = 0;
    _tt = -1;
    _aet = 0;
    
    int initialValue = sqrt(_lfsc).round();
    for (int i = 0; i < _BUFFER_SIZE; i++) {
      _ebuf[i] = initialValue;
      _lbuf[i] = 0;
    }
    
    // Установка порогов
    _t0 = T0;
    _ta = 3 * T0;
    _t1 = (_ta / 3).round();  // T1 = Ta/3 (вместо T0 как в оригинале после обучения)
    
    _isLearning = false;
    _lastDetection = -1000;
    _timer = 0;
    _inRefractory = false;
    _samplesProcessed = 0;
    _detections.clear();
    _jPoints.clear();
    _minutes = 0;
    _nextMinute = (60 * _sampleRate).round(); // 1 минута
    
    // Главный цикл обработки
    for (int t = 0; t < signal.length; t++) {
      bool detected = _processSample(signal[t], signal, t, false);
      
      // Эмуляция пропуска семплов во время рефрактерного периода
      if (detected) {
        // В оригинале: t += EyeClosing (пропуск семплов)
        // Здесь мы пропускаем итерации цикла
        int skipTo = t + _eyeClosingSamples;
        if (skipTo > t + 1) {
          // Пропускаем семплы, но нужно обновить состояние фильтра
          // Для этого вызываем _ltsamp для пропущенных семплов
          for (int tt = t + 1; tt < skipTo && tt < signal.length; tt++) {
            _ltsamp(tt, signal);
          }
          t = skipTo - 1; // -1 потому что цикл инкрементирует
          _inRefractory = false;
        }
      }
      
      // Отслеживание прогресса (как в оригинале - точки каждую минуту)
      if (t >= _nextMinute) {
        _nextMinute += (60 * _sampleRate).round();
        if (_debugMode) {
          stdout.write('.');
          if (++_minutes >= 60) {
            _minutes = 0;
          }
        }
      }
      
      _samplesProcessed++;
    }
    
    if (_debugMode) {
      if (_minutes > 0) print('');
      print('=== Detection Complete ===');
      print('Total detections: ${_detections.length}');
      
      double durationMinutes = signal.length / _sampleRate / 60.0;
      int expectedBeats = (durationMinutes * 75).round();
      print('Signal duration: ${durationMinutes.toStringAsFixed(1)} minutes');
      print('Expected beats: ~$expectedBeats (at 75 BPM)');
      
      if (_detections.length > expectedBeats * 2) {
        print('WARNING: Too many detections! Expected ~$expectedBeats, got ${_detections.length}');
        print('  Try increasing minThreshold');
      } else if (_detections.length < expectedBeats ~/ 2) {
        print('WARNING: Too few detections! Expected ~$expectedBeats, got ${_detections.length}');
        print('  Try decreasing minThreshold');
      } else {
        print('✓ Detection rate looks reasonable!');
      }
      
      if (_detections.isNotEmpty) {
        print('First 10 detections: ${_detections.take(10).toList()}');
        print('Last 10 detections: ${_detections.skip(max(0, _detections.length - 10)).toList()}');
      }
      print('==========================');
    }
    
    return List.from(_detections);
  }
  
  List<int> getJPoints() => List.from(_jPoints);
  int getThreshold() => _t1;
  bool isReady() => !_isLearning;
  int get detectionCount => _detections.length;
  
  Map<String, dynamic> getDebugInfo() {
    return {
      'detections': _detections.length,
      't0': _t0,
      'ta': _ta,
      't1': _t1,
      'isLearning': _isLearning,
      'lastDetection': _lastDetection,
      'gain': _gain,
      'inputIsRaw': _inputIsRaw,
    };
  }
}