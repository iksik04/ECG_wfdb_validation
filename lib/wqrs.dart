// lib/wqrs_detector.dart
import 'dart:math';

/// Класс для обнаружения QRS комплексов на основе алгоритма wqrs
/// Адаптирован для сигналов в милливольтах
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
  int _t0 = 0;
  int _ta = 0;
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
  
  // Параметры масштабирования
  final double _gain;  // Усиление в АЦП/мВ
  final bool _debugMode;
  int _samplesProcessed = 0;
  
  // Константы для преобразования
  static const int ADC_PER_MV = 200;  // Стандартное усиление ЭКГ
  
  WQRSDetector({
    required double sampleRate,
    double gain = 200.0,  // АЦП/мВ (по умолчанию 200)
    int powerLineFreq = 50,
    int minThreshold = 100,
    double eyeClosingPeriod = 0.25,
    double maxQRSWidth = 0.13,
    double ndp = 2.5,
    bool debugMode = false,
  }) : _sampleRate = sampleRate,
       _gain = gain,
       _minThreshold = minThreshold,
       _lpn = (sampleRate / powerLineFreq).floor().clamp(1, 8),
       _lp2n = (2 * (sampleRate / powerLineFreq).floor()).clamp(2, 16),
       _eyeClosingSamples = (eyeClosingPeriod * sampleRate).round(),
       _expectPeriod = (ndp * sampleRate).round(),
       _ltWindow = (maxQRSWidth * sampleRate).round().clamp(1, 256),
       _lfsc = 1.25 * gain * gain / sampleRate,
       _debugMode = debugMode {
    
    // Инициализация буфера
    int initialValue = sqrt(_lfsc).round();
    for (int i = 0; i < _BUFFER_SIZE; i++) {
      _ebuf[i] = initialValue;
    }
    
    if (_debugMode) {
      print('=== WQRS Detector Initialized ===');
      print('Sample Rate: $_sampleRate Hz');
      print('Gain: $_gain ADC/mV');
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
  }
  
  /// Преобразование сигнала из мВ в АЦП-единицы
  int _mvToAdc(double value) {
    return (value * _gain).round();
  }
  
  int _ltsamp(int t, List<double> signal) {
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
      
      // Получение отсчетов с преобразованием в АЦП
      v0 = (_tt >= 0 && _tt < signal.length) ? _mvToAdc(signal[_tt]) : 0;
      v1 = (_tt - _lpn >= 0 && _tt - _lpn < signal.length) ? _mvToAdc(signal[_tt - _lpn]) : 0;
      v2 = (_tt - _lp2n >= 0 && _tt - _lp2n < signal.length) ? _mvToAdc(signal[_tt - _lp2n]) : 0;
      
      // Фильтр низких частот
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
  
  bool processSample(double sample, List<double> signal) {
    int t = _samplesProcessed;
    
    int lt = _ltsamp(t, signal);
    _samplesProcessed++;
    
    // === ФАЗА ОБУЧЕНИЯ ===
    if (_isLearning) {
      _t0 += lt;
      
      int learningEnd = (8.0 * _sampleRate).round().clamp(1, _BUFFER_SIZE ~/ 2);
      
      if (t >= learningEnd) {
        _t0 ~/= (t + 1);
        _ta = 3 * _t0;
        _t1_learning = 2 * _t0;
        
        if (_debugMode) {
          print('=== Learning Complete ===');
          print('T0: $_t0');
          print('Ta: $_ta');
          print('T1 (learning): $_t1_learning');
          print('=========================');
        }
        
        _isLearning = false;
        return false;
      }
      
      return false;
    }
    
    // === ФАЗА ДЕТЕКЦИИ ===
    int currentT1 = _isLearning ? _t1_learning : (_ta ~/ 3);
    
    // Проверка рефрактерного периода
    if (_inRefractory) {
      if (t - _lastDetection >= _eyeClosingSamples) {
        _inRefractory = false;
      } else {
        return false;
      }
    }
    
    // Отладочный вывод
    if (_debugMode && t % 5000 == 0) {
      print('Sample $t: LT=$lt, T1=$currentT1, detections=${_detections.length}');
    }
    
    // Сравнение с порогом
    if (lt > currentT1) {
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
        // Поиск начала QRS
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
        
        // Регистрация обнаружения
        if (!_isLearning && tpq >= 0 && tpq < signal.length) {
          // Проверка, не слишком ли близко к предыдущему обнаружению
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
            _ta += (maxVal - _ta) ~/ 10;
            _t1 = _ta ~/ 3;
          }
        }
        
        return true;
      }
    } else {
      // Уменьшение порога
      _timer++;
      if (_timer > _expectPeriod && _ta > _minThreshold) {
        _ta--;
        _t1 = _ta ~/ 3;
        _timer = 0;
        
        if (_debugMode && _ta % 100 == 0) {
          print('Threshold lowered: Ta=$_ta, T1=$_t1');
        }
      }
    }
    
    return false;
  }
  
  List<int> detect(List<double> signal) {
    if (_debugMode) {
      print('=== Starting Detection ===');
      print('Signal length: ${signal.length} samples');
      
      if (signal.isNotEmpty) {
        double minSignal = signal.reduce((a, b) => a < b ? a : b);
        double maxSignal = signal.reduce((a, b) => a > b ? a : b);
        print('Signal range: ${minSignal.toStringAsFixed(2)} to ${maxSignal.toStringAsFixed(2)} mV');
        
        // Расчет ожидаемой амплитуды в АЦП
        double amplitude = maxSignal - minSignal;
        double expectedAdcAmplitude = amplitude * _gain;
        print('Expected ADC amplitude: ${expectedAdcAmplitude.toStringAsFixed(0)}');
        print('Recommended gain: ${(1000 / amplitude).toStringAsFixed(1)} ADC/mV');
      }
    }
    
    reset();
    
    // === ОБУЧЕНИЕ ===
    if (_debugMode) print('Phase 1: Learning...');
    
    int learningEnd = (8.0 * _sampleRate).round().clamp(1, _BUFFER_SIZE ~/ 2);
    int actualLearning = min(learningEnd, signal.length);
    
    for (int i = 0; i < actualLearning; i++) {
      _ltsamp(i, signal);
      int lt = _ltsamp(i, signal);
      _t0 += lt;
    }
    
    if (actualLearning > 0) {
      _t0 ~/= actualLearning;
      _ta = 3 * _t0;
      _t1_learning = 2 * _t0;
      _isLearning = false;
    } else {
      if (_debugMode) print('ERROR: Signal too short for learning');
      return [];
    }
    
    if (_debugMode) {
      print('Learning complete: T0=$_t0, Ta=$_ta, T1_learning=$_t1_learning');
    }
    
    // === ДЕТЕКЦИЯ ===
    if (_debugMode) print('Phase 2: Detection...');
    
    // Сброс состояния
    _yn = 0;
    _yn1 = 0;
    _yn2 = 0;
    _tt = -1;
    _aet = 0;
    _lastDetection = -1000;
    _timer = 0;
    _inRefractory = false;
    _samplesProcessed = 0;
    _detections.clear();
    
    int initialValue = sqrt(_lfsc).round();
    for (int i = 0; i < _BUFFER_SIZE; i++) {
      _ebuf[i] = initialValue;
      _lbuf[i] = 0;
    }
    
    // Обработка сигнала
    for (int i = 0; i < signal.length; i++) {
      processSample(signal[i], signal);
      
      if (_debugMode && i % (signal.length ~/ 10) == 0 && i > 0) {
        print('Progress: ${(i / signal.length * 100).toStringAsFixed(0)}% '
              '(${_detections.length} detections)');
      }
    }
    
    if (_debugMode) {
      print('=== Detection Complete ===');
      print('Total detections: ${_detections.length}');
      
      // Расчет ожидаемого количества QRS
      double durationMinutes = signal.length / _sampleRate / 60.0;
      int expectedBeats = (durationMinutes * 75).round(); // 75 ударов в минуту
      print('Signal duration: ${durationMinutes.toStringAsFixed(1)} minutes');
      print('Expected beats: ~$expectedBeats (at 75 BPM)');
      
      if (_detections.length > expectedBeats * 2) {
        print('WARNING: Too many detections! Expected ~$expectedBeats, got ${_detections.length}');
        print('  Try increasing gain or minThreshold');
      } else if (_detections.length < expectedBeats ~/ 2) {
        print('WARNING: Too few detections! Expected ~$expectedBeats, got ${_detections.length}');
        print('  Try decreasing gain or minThreshold');
      } else {
        print('Detection rate looks reasonable!');
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
    };
  }
}