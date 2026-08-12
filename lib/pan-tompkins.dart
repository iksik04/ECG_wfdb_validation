// lib/pan_tompkins_alg.dart
import 'dart:math';

// ============================================================
// Класс PanTompkinsQRS (основная логика детекции и обработки)
// ============================================================
class PanTompkinsQRS {
  // Буферы для ФНЧ
  final List<double> _lBuffer = [];
  double _lPrev1 = 0.0;
  double _lPrev2 = 0.0;

  // Буферы для ФВЧ
  final List<double> _hBuffer = [];
  double _hPrev1 = 0.0;

  double _sampleRate = 250.0;

  // Состояния для обучения
  bool _isLearningPhase1 = false;
  bool _isLearningPhase2 = false;
  bool _isInitialized = false;
  int _learningPhase1Counter = 0;
  int _learningPhase2Counter = 0;
  List<double> _learningPeaks = [];
  List<int> _learningPeakIndices = [];
  
  // Время последнего обнаруженного QRS (в индексах) для рефрактерного периода
  int _lastQRSIndex = -1000;
  
  // Предыдущий наклон для T-волновой дискриминации
  double _prevSlope = 0.0;

  // Адаптивные пороги для интегрального сигнала
  double _SPKI = 0.0;
  double _NPKI = 0.0;
  double _THRESHOLDI1 = 0.0;
  double _THRESHOLDI2 = 0.0;

  // Адаптивные пороги для фильтрованного сигнала
  double _SPKF = 0.0;
  double _NPKF = 0.0;
  double _THRESHOLDF1 = 0.0;
  double _THRESHOLDF2 = 0.0;

  // Оценки RR-интервалов
  final List<double> _rrAvg1 = [];
  final List<double> _rrAvg2 = [];

  double _RRLowLimit = 0.0;
  double _RRHighLimit = 0.0;
  double _RRMissedLimit = 0.0;

  // Флаг для отслеживания, был ли уже обновлен RR_avg2 в текущем цикле
  bool _rrAvg2Updated = false;

  PanTompkinsQRS({double sampleRate = 250.0}) {
    _sampleRate = sampleRate;
    reset();
  }

  void reset() {
    _lBuffer.clear();
    _lPrev1 = 0.0;
    _lPrev2 = 0.0;
    _hBuffer.clear();
    _hPrev1 = 0.0;
    
    _isLearningPhase1 = false;
    _isLearningPhase2 = false;
    _isInitialized = false;
    _learningPhase1Counter = 0;
    _learningPhase2Counter = 0;
    _learningPeaks.clear();
    _learningPeakIndices.clear();
    
    _lastQRSIndex = -1000;
    _prevSlope = 0.0;
    
    _SPKI = 0.0;
    _NPKI = 0.0;
    _THRESHOLDI1 = 0.0;
    _THRESHOLDI2 = 0.0;
    _SPKF = 0.0;
    _NPKF = 0.0;
    _THRESHOLDF1 = 0.0;
    _THRESHOLDF2 = 0.0;
    _rrAvg1.clear();
    _rrAvg2.clear();
    _RRLowLimit = 0.0;
    _RRHighLimit = 0.0;
    _RRMissedLimit = 0.0;
    _rrAvg2Updated = false;
  }

  /// ФНЧ: y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[n-8] + x[n-16]
  double _applyLowPassFilter_250(double val) {
    _lBuffer.add(val);
    if (_lBuffer.length > 17) _lBuffer.removeAt(0);

    double xn = _lBuffer[_lBuffer.length - 1];
    double xn8 = (_lBuffer.length - 9 >= 0) ? _lBuffer[_lBuffer.length - 9] : 0.0;
    double xn16 = (_lBuffer.length - 17 >= 0) ? _lBuffer[_lBuffer.length - 17] : 0.0;

    double y = 2.0 * _lPrev1 - _lPrev2 + xn - 2.0 * xn8 + xn16;

    _lPrev2 = _lPrev1;
    _lPrev1 = y;
    return y;
  }

  /// ФВЧ: y[n] = 40*x[n-20] - y[n-1] - x[n] + x[n-40]
  double _applyHighPassFilter_250(double val) {
    _hBuffer.add(val);
    if (_hBuffer.length > 41) _hBuffer.removeAt(0);

    double xn = _hBuffer[_hBuffer.length - 1];
    double xn20 = (_hBuffer.length - 21 >= 0) ? _hBuffer[_hBuffer.length - 21] : 0.0;
    double xn40 = (_hBuffer.length - 41 >= 0) ? _hBuffer[_hBuffer.length - 41] : 0.0;

    double y = 40.0 * xn20 - _hPrev1 - xn + xn40;

    _hPrev1 = y;
    return y;
  }

  /// ФНЧ: y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[n-11] + x[n-22]
  double _applyLowPassFilter_360(double val) {
    _lBuffer.add(val);
    if (_lBuffer.length > 23) _lBuffer.removeAt(0);

    double xn = _lBuffer[_lBuffer.length - 1];
    double xn11 = (_lBuffer.length - 12 >= 0) ? _lBuffer[_lBuffer.length - 12] : 0.0;
    double xn22 = (_lBuffer.length - 23 >= 0) ? _lBuffer[_lBuffer.length - 23] : 0.0;

    double y = 2.0 * _lPrev1 - _lPrev2 + xn - 2.0 * xn11 + xn22;

    _lPrev2 = _lPrev1;
    _lPrev1 = y;
    return y;
  }

  /// ФВЧ: y[n] = 58*x[n-29] - y[n-1] - x[n] + x[n-58]
  double _applyHighPassFilter_360(double val) {
    _hBuffer.add(val);
    if (_hBuffer.length > 59) _hBuffer.removeAt(0);

    double xn = _hBuffer[_hBuffer.length - 1];
    double xn29 = (_hBuffer.length - 30 >= 0) ? _hBuffer[_hBuffer.length - 30] : 0.0;
    double xn58 = (_hBuffer.length - 59 >= 0) ? _hBuffer[_hBuffer.length - 59] : 0.0;

    double y = 58.0 * xn29 - _hPrev1 - xn + xn58;

    _hPrev1 = y;
    return y;
  }

  // Дифференцирование
  List<double> derivative(List<double> signal, double fs) {
    final int n = signal.length;
    final List<double> result = List<double>.filled(n, 0.0);
    final double scale = fs / 8.0;

    if (n > 1) {
      result[1] = -2.0 * signal[0] * scale;
    }

    for (int i = 2; i <= n - 3; i++) {
      double val = -2.0 * signal[i - 1] - signal[i - 2] + 2.0 * signal[i + 1] + signal[i + 2];
      result[i] = val * scale;
    }

    if (n >= 3) {
      int i = n - 2;
      result[i] = (-2.0 * signal[i - 1] - signal[i - 2]) * scale;
    }
    if (n >= 4) {
      int i = n - 1;
      result[i] = (-2.0 * signal[i - 1] - signal[i - 2]) * scale;
    }
    return result;
  }

  List<double> squaring(List<double> signal) {
    final int n = signal.length;
    final List<double> result = List<double>.filled(n, 0.0);

    for (int index = 0; index < signal.length; index++) {
      result[index] = signal[index] * signal[index];
    }

    return result;
  }

  List<double> movingWindowIntegration(List<double> signal, double fs) {
    final int n = signal.length;
    final int winSize = (0.150 * fs).round();
    final int effectiveWin = winSize > n ? n : winSize;
    final double invWin = 1.0 / winSize;

    final List<double> result = List<double>.filled(n, 0.0);
    double sumRaw = 0.0;

    for (int j = 0; j < effectiveWin; j++) {
      sumRaw += signal[j];
      result[j] = sumRaw * invWin;
    }

    for (int i = winSize; i < n; i++) {
      sumRaw += signal[i] - signal[i - winSize];
      result[i] = sumRaw * invWin;
    }
    return result;
  }

  // ---- Вспомогательные методы ----

  List<int> _findLocalMaxima(List<double> signal) {
    final List<int> peaks = [];
    for (int i = 1; i < signal.length - 1; i++) {
      if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1]) {
        peaks.add(i);
      }
    }
    return peaks;
  }

  void _updateThresholds() {
    _THRESHOLDI1 = _NPKI + 0.25 * (_SPKI - _NPKI);
    _THRESHOLDI2 = 0.5 * _THRESHOLDI1;
    _THRESHOLDF1 = _NPKF + 0.25 * (_SPKF - _NPKF);
    _THRESHOLDF2 = 0.5 * _THRESHOLDF1;
  }

  /// ИСПРАВЛЕННАЯ ЛОГИКА ОБНОВЛЕНИЯ RR-ИНТЕРВАЛОВ
  /// В соответствии со статьей Pan & Tompkins:
  /// 
  /// RR_avg1 = среднее 8 последних RR-интервалов (любых)
  /// RR_avg2 = среднее 8 последних RR-интервалов, которые попадают в допустимый диапазон (92-116% от текущего RR_avg2)
  /// 
  /// Если все 8 последних интервалов из RR_avg1 попадают в допустимый диапазон,
  /// то RR_avg2 = RR_avg1 (обновление RR_avg2)
  void _updateRRIntervals(double rr) {
    // 1. Обновляем RR_avg1 (всегда добавляем новый интервал)
    _rrAvg1.add(rr);
    if (_rrAvg1.length > 8) _rrAvg1.removeAt(0);

    // 2. Обновляем RR_avg2 только если RR_avg2 уже инициализирован
    if (_rrAvg2.isNotEmpty) {
      // Проверяем, попадает ли новый RR-интервал в допустимый диапазон
      if (rr >= _RRLowLimit && rr <= _RRHighLimit) {
        // Добавляем интервал в RR_avg2
        _rrAvg2.add(rr);
        if (_rrAvg2.length > 8) _rrAvg2.removeAt(0);
      }
      
      // 3. Проверяем, все ли 8 последних интервалов из RR_avg1 попадают в допустимый диапазон
      // Это соответствует формуле (29) из статьи: 
      // "If each of the eight most-recent sequential RR intervals that are calculated 
      // from RR AVERAGE1 is between the RR LOW LIMIT and the RR HIGH LIMIT, 
      // we interpret the heart rate to be regular for these eight heart beats and 
      // RR AVERAGE2 <- RR AVERAGE1."
      if (_rrAvg1.length >= 8) {
        bool allInRange = true;
        for (int i = 0; i < _rrAvg1.length; i++) {
          if (_rrAvg1[i] < _RRLowLimit || _rrAvg1[i] > _RRHighLimit) {
            allInRange = false;
            break;
          }
        }
        
        // Если все 8 интервалов в допустимом диапазоне, обновляем RR_avg2 = RR_avg1
        if (allInRange) {
          final double avg1 = _rrAvg1.reduce((a, b) => a + b) / _rrAvg1.length;
          // Очищаем RR_avg2 и заполняем его копией RR_avg1
          _rrAvg2.clear();
          _rrAvg2.addAll(_rrAvg1);
          // Убеждаемся, что длина не превышает 8
          while (_rrAvg2.length > 8) {
            _rrAvg2.removeAt(0);
          }
          _rrAvg2Updated = true;
        }
      }
    } else {
      // Если RR_avg2 пуст, добавляем интервал (первая инициализация)
      _rrAvg2.add(rr);
    }

    // 4. Пересчитываем лимиты на основе текущего RR_avg2
    if (_rrAvg2.length >= 8) {
      final double avg2 = _rrAvg2.reduce((a, b) => a + b) / _rrAvg2.length;
      _RRLowLimit = 0.92 * avg2;
      _RRHighLimit = 1.16 * avg2;
      _RRMissedLimit = 1.66 * avg2;
    } else if (_rrAvg2.isNotEmpty) {
      // Если меньше 8 интервалов, используем доступные для приблизительной оценки
      final double avg2 = _rrAvg2.reduce((a, b) => a + b) / _rrAvg2.length;
      _RRLowLimit = 0.92 * avg2;
      _RRHighLimit = 1.16 * avg2;
      _RRMissedLimit = 1.66 * avg2;
    }
  }

  int _performSearchBack({
    required List<double> integrationSignal,
    required List<double> bandPassSignal,
    required int currentIntegIdx,
    required int fs,
  }) {
    final int searchWindow = (_RRMissedLimit * fs).round();
    int start = (currentIntegIdx - searchWindow).clamp(0, integrationSignal.length - 1);
    int end = currentIntegIdx.clamp(0, integrationSignal.length - 1);

    int bestIntegIdx = -1;
    double maxVal = -double.infinity;
    for (int i = start; i <= end; i++) {
      final double val = integrationSignal[i];
      if (val > _THRESHOLDI2 && val < _THRESHOLDI1 && val > maxVal) {
        maxVal = val;
        bestIntegIdx = i;
      }
    }
    if (bestIntegIdx == -1) return -1;

    int leftF = (bestIntegIdx - (0.15 * fs).round()).clamp(0, bandPassSignal.length - 1);
    int rightF = (bestIntegIdx + (0.15 * fs).round() + 1).clamp(0, bandPassSignal.length);
    int bestFiltIdx = -1;
    double maxFilt = -double.infinity;
    for (int j = leftF; j < rightF; j++) {
      if (bandPassSignal[j] > _THRESHOLDF2 && bandPassSignal[j] > maxFilt) {
        maxFilt = bandPassSignal[j];
        bestFiltIdx = j;
      }
    }
    if (bestFiltIdx == -1) return -1;

    return bestIntegIdx;
  }

  int _findFiltPeak(List<double> bandPassSignal, int integIdx, int fs) {
    final int window = (0.15 * fs).round();
    int left = (integIdx - window).clamp(0, bandPassSignal.length - 1);
    int right = (integIdx + window + 1).clamp(0, bandPassSignal.length);
    double maxVal = -double.infinity;
    int maxIdx = -1;
    for (int j = left; j < right; j++) {
      if (bandPassSignal[j] > maxVal) {
        maxVal = bandPassSignal[j];
        maxIdx = j;
      }
    }
    return maxIdx;
  }

  double _computeMaxSlope(List<double> signal, int idx, int fs) {
    final int half = (0.075 * fs).round();
    int start = (idx - half).clamp(0, signal.length - 1);
    int end = (idx + half).clamp(0, signal.length - 1);
    double maxSlope = 0.0;
    for (int i = start; i < end; i++) {
      final double diff = signal[i + 1] - signal[i];
      if (diff > maxSlope) maxSlope = diff;
    }
    return maxSlope;
  }

  /// Масштабирование обработанного сигнала к уровню исходного
  List<double> scaleToOriginal(List<double> processed, List<double> original) {
    if (processed.length != original.length) {
      throw ArgumentError('Длины сигналов должны совпадать.');
    }

    double rmsOrig = 0.0;
    for (var v in original) rmsOrig += v * v;
    rmsOrig = sqrt(rmsOrig / original.length);

    double rmsProc = 0.0;
    for (var v in processed) rmsProc += v * v;
    rmsProc = sqrt(rmsProc / processed.length);

    if (rmsProc == 0.0) return List.from(processed);

    double scale = rmsOrig / rmsProc;
    return processed.map((v) => v * scale).toList();
  }

  /// Вычисление задержки в зависимости от частоты дискретизации
  int _getDelay(int fs) {
    if (fs == 250) {
      return 24;
    } else if (fs == 360) {
      return 42;
    } else return 20;
  }

  // ============================================================
  // ДВУХФАЗНОЕ ОБУЧЕНИЕ
  // ============================================================

  /// Основной метод обработки сигнала с двухфазным обучением
  List<double> process(List<double> signal) {
    if (signal.isEmpty) return [];

    // Сброс всех внутренних состояний
    reset();

    // 1. Полосовая фильтрация
    final List<double> filtered = List.filled(signal.length, 0.0);
    if (_sampleRate == 250) {
      for (int i = 0; i < signal.length; i++) {
        double low = _applyLowPassFilter_250(signal[i]);
        filtered[i] = _applyHighPassFilter_250(low) / 2560;
      }
    } else if (_sampleRate == 360) {
      for (int i = 0; i < signal.length; i++) {
        double low = _applyLowPassFilter_360(signal[i]);
        filtered[i] = _applyHighPassFilter_360(low) / 7018;
      }
    }

    // 2. Производная
    final List<double> derivated = derivative(filtered, _sampleRate);

    // 3. Возведение в квадрат
    final List<double> squared = squaring(derivated);

    // 4. Скользящее интегрирование
    final List<double> integrated = movingWindowIntegration(squared, _sampleRate);

    // 5. Двухфазное обучение
    _performTwoPhaseLearning(integrated, filtered);

    // 6. Детекция пиков с использованием обученных порогов
    final List<int> peaks = detectPeaks(
      ecgSignal: signal,
      fs: _sampleRate.round(),
      integrationSignal: integrated,
      bandPassSignal: filtered,
    );

    // 7. Масштабирование
    final List<double> scaled = scaleToOriginal(integrated, signal);

    // 8. Коррекция задержки
    int delay = _getDelay(_sampleRate.round());
    final List<double> shifted = List.filled(scaled.length, 0.0);
    final int copyLength = scaled.length - delay;
    for (int i = 0; i < copyLength; i++) {
      shifted[i] = scaled[i + delay];
    }

    return shifted;
  }

  /// Двухфазное обучение по статье Pan & Tompkins
  void _performTwoPhaseLearning(
    List<double> integrationSignal,
    List<double> bandPassSignal,
  ) {
    if (integrationSignal.length < 100) return;

    // Находим все локальные максимумы в интегральном сигнале
    final List<int> integPeaks = _findLocalMaxima(integrationSignal);
    if (integPeaks.length < 2) return;

    // === ФАЗА 1: Обучение (первые ~2 секунды) ===
    _isLearningPhase1 = true;
    _learningPhase1Counter = 0;
    _learningPeaks.clear();
    _learningPeakIndices.clear();

    // Используем первые 2 секунды сигнала для инициализации порогов
    final int learningWindowSamples = (2.0 * _sampleRate).round();
    int learningEnd = min(learningWindowSamples, integPeaks.length);

    // Находим максимальные пики в первом окне обучения
    double maxIntegPeak = 0.0;
    double maxFiltPeak = 0.0;
    int maxIntegIdx = -1;
    int maxFiltIdx = -1;

    for (int i = 0; i < learningEnd && i < integPeaks.length; i++) {
      int idx = integPeaks[i];
      if (idx < integrationSignal.length) {
        if (integrationSignal[idx] > maxIntegPeak) {
          maxIntegPeak = integrationSignal[idx];
          maxIntegIdx = idx;
          
          // Находим соответствующий пик в фильтрованном сигнале
          int filtIdx = _findFiltPeak(bandPassSignal, idx, _sampleRate.round());
          if (filtIdx != -1 && filtIdx < bandPassSignal.length) {
            if (bandPassSignal[filtIdx] > maxFiltPeak) {
              maxFiltPeak = bandPassSignal[filtIdx];
              maxFiltIdx = filtIdx;
            }
          }
        }
      }
    }

    // Инициализация порогов на основе максимального пика
    if (maxIntegIdx != -1 && maxFiltIdx != -1) {
      _SPKI = maxIntegPeak;
      _SPKF = maxFiltPeak;
      _NPKI = 0.5 * maxIntegPeak;
      _NPKF = 0.5 * maxFiltPeak;
      _updateThresholds();
      
      // Запоминаем первый обнаруженный QRS
      _lastQRSIndex = maxIntegIdx;
      _prevSlope = _computeMaxSlope(bandPassSignal, maxFiltIdx, _sampleRate.round());
      
      // Добавляем в список обученных пиков
      _learningPeaks.add(maxIntegPeak);
      _learningPeakIndices.add(maxIntegIdx);
    }

    _isLearningPhase1 = false;
    _isLearningPhase2 = true;

    // === ФАЗА 2: Обучение на двух последовательных QRS ===
    _learningPhase2Counter = 0;
    int rr1 = -1;
    int rr2 = -1;

    // Находим два последовательных QRS комплекса для инициализации RR-интервалов
    for (int i = 0; i < integPeaks.length && _learningPhase2Counter < 2; i++) {
      int idx = integPeaks[i];
      if (idx < integrationSignal.length) {
        // Проверяем, что это валидный QRS
        double integVal = integrationSignal[idx];
        int filtIdx = _findFiltPeak(bandPassSignal, idx, _sampleRate.round());
        double filtVal = filtIdx != -1 ? bandPassSignal[filtIdx] : 0.0;
        
        if (integVal >= _THRESHOLDI1 && filtVal >= _THRESHOLDF1) {
          if (_learningPhase2Counter == 0) {
            rr1 = idx;
            _learningPhase2Counter++;
          } else if (_learningPhase2Counter == 1) {
            rr2 = idx;
            // Вычисляем RR-интервал
            double rrInterval = (rr2 - rr1).abs() / _sampleRate;
            if (rrInterval > 0.2 && rrInterval < 2.0) { // Физиологический диапазон
              _updateRRIntervals(rrInterval);
              _lastQRSIndex = idx;
              _prevSlope = _computeMaxSlope(bandPassSignal, 
                _findFiltPeak(bandPassSignal, idx, _sampleRate.round()), 
                _sampleRate.round());
            }
            _learningPhase2Counter++;
          }
        }
      }
    }

    _isLearningPhase2 = false;
    _isInitialized = true;
  }

  /// Основной метод детекции с поддержкой обучения
  List<int> detectPeaks({
    required List<double> ecgSignal,
    required int fs,
    required List<double> integrationSignal,
    required List<double> bandPassSignal,
  }) {
    final int window = (0.15 * fs).round();
    final List<int> integPeaks = _findLocalMaxima(integrationSignal);

    final List<_PeakPair> candidatePairs = [];

    for (int idx in integPeaks) {
      int left = (idx - window).clamp(0, bandPassSignal.length - 1);
      int right = (idx + window + 1).clamp(0, bandPassSignal.length);
      double maxVal = -double.infinity;
      int maxIdx = -1;
      for (int j = left; j < right; j++) {
        if (bandPassSignal[j] > maxVal) {
          maxVal = bandPassSignal[j];
          maxIdx = j;
        }
      }
      if (maxIdx != -1) {
        candidatePairs.add(_PeakPair(integIdx: idx, filtIdx: maxIdx));
      }
    }

    final List<int> rPeaks = [];

    // Если обучение не завершено, пропускаем детекцию
    if (!_isInitialized) {
      return rPeaks;
    }

    for (int i = 0; i < candidatePairs.length; i++) {
      final pair = candidatePairs[i];
      final int integIdx = pair.integIdx;
      final int filtIdx = pair.filtIdx;
      final double integVal = integrationSignal[integIdx];
      final double filtVal = bandPassSignal[filtIdx];

      // Рефрактерный период (200 мс)
      if (integIdx - _lastQRSIndex < (0.2 * fs).round()) {
        continue;
      }

      // Обновление RR-интервалов (только если есть предыдущий QRS)
      if (i > 0 && _lastQRSIndex > 0) {
        final double rr = (integIdx - _lastQRSIndex) / fs;
        if (rr > 0.2 && rr < 2.0) { // Физиологический диапазон
          _updateRRIntervals(rr);
        }
      }

      // Проверка на нерегулярный ритм
      bool irregularRate = false;
      if (_rrAvg1.length >= 8) {
        irregularRate = true;
        for (int j = _rrAvg1.length - 8; j < _rrAvg1.length; j++) {
          if (_rrAvg1[j] >= _RRLowLimit && _rrAvg1[j] <= _RRHighLimit) {
            irregularRate = false;
            break;
          }
        }
      }

      double thresholdI1 = irregularRate ? _THRESHOLDI1 * 0.5 : _THRESHOLDI1;
      double thresholdF1 = irregularRate ? _THRESHOLDF1 * 0.5 : _THRESHOLDF1;

      bool isQRS = false;

      // Основная проверка
      if (integVal >= thresholdI1 && filtVal >= thresholdF1) {
        isQRS = true;
        _SPKI = 0.125 * integVal + 0.875 * _SPKI;
        _SPKF = 0.125 * filtVal + 0.875 * _SPKF;
      } 
      // Search-back для пропущенных QRS
      else if (_rrAvg2.isNotEmpty &&
          i > 0 &&
          (integIdx - _lastQRSIndex) / fs > _RRMissedLimit) {
        final int searchBack = _performSearchBack(
          integrationSignal: integrationSignal,
          bandPassSignal: bandPassSignal,
          currentIntegIdx: integIdx,
          fs: fs,
        );
        if (searchBack != -1) {
          isQRS = true;
          int filtSearchIdx = _findFiltPeak(bandPassSignal, searchBack, fs);
          if (filtSearchIdx != -1) {
            _SPKI = 0.25 * integrationSignal[searchBack] + 0.75 * _SPKI;
            _SPKF = 0.25 * bandPassSignal[filtSearchIdx] + 0.75 * _SPKF;
            rPeaks.add(searchBack);
            _lastQRSIndex = searchBack;
            _updateThresholds();
            continue;
          }
        }
      }

      if (isQRS) {
        // T-волновая дискриминация
        final double rr = i > 0 ? (integIdx - _lastQRSIndex) / fs : 0;
        if (rr > 0.2 && rr < 0.36) {
          final double currSlope = _computeMaxSlope(bandPassSignal, filtIdx, fs);
          if (currSlope < 0.5 * _prevSlope) {
            _NPKI = 0.125 * integVal + 0.875 * _NPKI;
            _NPKF = 0.125 * filtVal + 0.875 * _NPKF;
            _updateThresholds();
            continue;
          }
        }

        rPeaks.add(integIdx);
        _lastQRSIndex = integIdx;
        _prevSlope = _computeMaxSlope(bandPassSignal, filtIdx, fs);
        _updateThresholds();
      } else {
        _NPKI = 0.125 * integVal + 0.875 * _NPKI;
        _NPKF = 0.125 * filtVal + 0.875 * _NPKF;
        _updateThresholds();
      }
    }

    return rPeaks;
  }

  /// Дополнительный метод для получения R-пиков
  List<int> detectRPeaks(List<double> signal) {
    if (signal.isEmpty) return [];

    reset();

    final List<double> filtered = List.filled(signal.length, 0.0);
    if (_sampleRate == 250) {
      for (int i = 0; i < signal.length; i++) {
        double low = _applyLowPassFilter_250(signal[i]);
        filtered[i] = _applyHighPassFilter_250(low) / 2560;
      }
    } else if (_sampleRate == 360) {
      for (int i = 0; i < signal.length; i++) {
        double low = _applyLowPassFilter_360(signal[i]);
        filtered[i] = _applyHighPassFilter_360(low) / 7018;
      }
    }

    final List<double> derivated = derivative(filtered, _sampleRate);
    final List<double> squared = squaring(derivated);
    final List<double> integrated = movingWindowIntegration(squared, _sampleRate);

    // Выполняем обучение
    _performTwoPhaseLearning(integrated, filtered);

    final List<int> peaks = detectPeaks(
      ecgSignal: signal,
      fs: _sampleRate.round(),
      integrationSignal: integrated,
      bandPassSignal: filtered,
    );

    int delay = _getDelay(_sampleRate.round());
    return peaks.map((p) => p - delay).where((p) => p >= 0).toList();
  }
}

// Вспомогательный класс для хранения пары индексов
class _PeakPair {
  final int integIdx;
  final int filtIdx;
  _PeakPair({required this.integIdx, required this.filtIdx});
}