// lib/pan_tompkins_alg.dart
import 'dart:math';

class PanTompkinsQRS {

  List<double> scaleToOriginal(List<double> processed, List<double> original) {
    if (processed.length != original.length) {
      throw ArgumentError('Длины сигналов должны совпадать.');
    }

    // Вычисляем RMS исходного сигнала
    double rmsOrig = 0.0;
    for (var v in original) rmsOrig += v * v;
    rmsOrig = sqrt(rmsOrig / original.length);

    // Вычисляем RMS обработанного сигнала
    double rmsProc = 0.0;
    for (var v in processed) rmsProc += v * v;
    rmsProc = sqrt(rmsProc / processed.length);

    // Если RMS обработанного сигнала равен нулю, возвращаем копию без изменений
    if (rmsProc == 0.0) return List.from(processed);

    double scale = rmsOrig / rmsProc;
    return processed.map((v) => v * scale).toList();
  }

  // Буферы для ФНЧ (нужно 13 отсчётов для x[n-12])
  final List<double> _lBuffer = [];
  double _lPrev1 = 0.0;  // y[n-1]
  double _lPrev2 = 0.0;  // y[n-2]

  // Буферы для ФВЧ (нужно 33 отсчёта для x[n-32])
  final List<double> _hBuffer = [];
  double _hPrev1 = 0.0;  // y[n-1]

  double _sampleRate = 250.0;

  PanTompkinsQRS({double sampleRate = 250.0}) {
    _sampleRate = sampleRate;
  }

  void reset() {
    _lBuffer.clear();
    _lPrev1 = 0.0;
    _lPrev2 = 0.0;
    _hBuffer.clear();
    _hPrev1 = 0.0;
  }

  /// ФНЧ: y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[n-6] + x[n-12]
  double _applyLowPassFilter(double val) {
    _lBuffer.add(val);
    if (_lBuffer.length > 13) _lBuffer.removeAt(0);

    double xn = _lBuffer[_lBuffer.length - 1];
    double xn6 = (_lBuffer.length - 7 >= 0) ? _lBuffer[_lBuffer.length - 7] : 0.0;
    double xn12 = (_lBuffer.length - 13 >= 0) ? _lBuffer[_lBuffer.length - 13] : 0.0;

    double y = 2.0 * _lPrev1 - _lPrev2 + xn - 2.0 * xn6 + xn12;

    _lPrev2 = _lPrev1;
    _lPrev1 = y;
    return y;
  }

  /// ФВЧ: y[n] = 32*x[n-16] - y[n-1] - x[n] + x[n-32]
  double _applyHighPassFilter(double val) {
    _hBuffer.add(val);
    if (_hBuffer.length > 33) _hBuffer.removeAt(0);

    double xn = _hBuffer[_hBuffer.length - 1];
    double xn16 = (_hBuffer.length - 17 >= 0) ? _hBuffer[_hBuffer.length - 17] : 0.0;
    double xn32 = (_hBuffer.length - 33 >= 0) ? _hBuffer[_hBuffer.length - 33] : 0.0;

    double y = 32.0 * xn16 - _hPrev1 - xn + xn32;

    _hPrev1 = y;
    return y;
  }

  // Дифференцирование
  List<double> derivative(List<double> signal, double fs) {
      final int n = signal.length;
      final List<double> result = List<double>.filled(n, 0.0);
      final double scale = fs / 8.0;

      // Индексы 0 и 1 – особые случаи (нет полного набора соседей)
      // index = 0: результат всегда 0 (оставляем как есть)

      if (n > 1) {
        // index = 1: только -2*signal[0]
        result[1] = -2.0 * signal[0] * scale;
      }

      // Основная часть: от 2 до n-3 включительно (если n >= 5)
      for (int i = 2; i <= n - 3; i++) {
        double val = -2.0 * signal[i - 1] - signal[i - 2] + 2.0 * signal[i + 1] + signal[i + 2];
        result[i] = val * scale;
      }

      // Индексы n-2 и n-1 (если они не вошли в основную часть)
      if (n >= 3) {
        // index = n-2: условия: index>=1, index>=2, index<=n-2 (да), index<=n-3? Нет, если n-2 > n-3, т.е. всегда нет.
        // Поэтому только отрицательные члены: -2*signal[i-1] - signal[i-2]
        int i = n - 2;
        result[i] = (-2.0 * signal[i - 1] - signal[i - 2]) * scale;
      }
      if (n >= 4) {
        // index = n-1: только условия index>=1 и index>=2 (отрицательные члены), положительных нет (т.к. i+1 >= n)
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
    print(fs);
    final int winSize = (0.150 * fs).round();
    // Если окно больше сигнала – обработаем все элементы как накопленное среднее
    final int effectiveWin = winSize > n ? n : winSize;
    final double invWin = 1.0 / winSize; // или 1.0 / effectiveWin, если менять семантику

    final List<double> result = List<double>.filled(n, 0.0);
    double sumRaw = 0.0;

    // Первые элементы – накопленное среднее (как в исходном коде)
    for (int j = 0; j < effectiveWin; j++) {
      sumRaw += signal[j];
      result[j] = sumRaw * invWin;
    }

    // Основной цикл скользящего окна (только если сигнал длиннее окна)
    for (int i = winSize; i < n; i++) {
      sumRaw += signal[i] - signal[i - winSize];
      result[i] = sumRaw * invWin;
    }
    return result;
  }

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
    final List<double> _rrAvg1 = []; // последние 8 RR (в секундах)
    final List<double> _rrAvg2 = []; // последние 8 RR в допустимом диапазоне

    double _RRLowLimit = 0.0;
    double _RRHighLimit = 0.0;
    double _RRMissedLimit = 0.0;

    // Предыдущий наклон для T-волновой дискриминации
    double _prevSlope = 0.0;

    // Время последнего обнаруженного QRS (в индексах) для рефрактерного периода
    int _lastQRSIndex = -1000;

    // Сброс всех порогов в начальное состояние (например, при смене пациента)
    void resetThresholds() {
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
      _prevSlope = 0.0;
      _lastQRSIndex = -1000;
    }

    // Основной метод детекции
    List<int> detectPeaks({
      required List<double> ecgSignal,
      required int fs,
      required List<double> integrationSignal,
      required List<double> bandPassSignal,
    }) {
      // 1. Находим все локальные максимумы в интегральном сигнале
      final int window = (0.15 * fs).round(); // ширина окна для поиска пика в фильтрованном
      final List<int> integPeaks = _findLocalMaxima(integrationSignal);

      // Список для хранения пар (индекс в интегральном, индекс в фильтрованном)
      final List<_PeakPair> candidatePairs = [];

      // Для каждого пика интегрального сигнала ищем соответствующий пик в фильтрованном
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

      // Список финальных fiducial меток (индексы в интегральном сигнале)
      final List<int> rPeaks = [];

      // Проходим по кандидатам
      for (int i = 0; i < candidatePairs.length; i++) {
        final pair = candidatePairs[i];
        final int integIdx = pair.integIdx;
        final int filtIdx = pair.filtIdx;

        // ---- Инициализация (первые два комплекса) ----
        if (i < 2) {
          // Первые два пика считаем QRS (обучение)
          _SPKI = 0.5 * integrationSignal[integIdx];
          _NPKI = 0.5 * _SPKI;
          _SPKF = 0.5 * bandPassSignal[filtIdx];
          _NPKF = 0.5 * _SPKF;
          _updateThresholds();
          rPeaks.add(integIdx);
          _lastQRSIndex = integIdx;
          // Запоминаем наклон для T-волн
          _prevSlope = _computeMaxSlope(integrationSignal, integIdx, fs);
          continue;
        }

        // ---- Проверка рефрактерного периода (200 мс) ----
        if (integIdx - _lastQRSIndex < (0.2 * fs).round()) {
          continue; // пропускаем пик
        }

        // ---- Обновление RR-интервалов и пределов ----
        _updateRRIntervals(candidatePairs, i, fs);

        // ---- Адаптация порогов при нерегулярном ритме ----
        if (_rrAvg1.isNotEmpty && _rrAvg1.last < _RRLowLimit) {
          _THRESHOLDI1 *= 0.5;
          _THRESHOLDF1 *= 0.5;
          _THRESHOLDI2 = 0.5 * _THRESHOLDI1;
          _THRESHOLDF2 = 0.5 * _THRESHOLDF1;
        }

        // ---- Проверка превышения порогов ----
        final double integVal = integrationSignal[integIdx];
        final double filtVal = bandPassSignal[filtIdx];

        bool isQRS = false;

        // Основной порог (первый)
        if (integVal >= _THRESHOLDI1 && filtVal >= _THRESHOLDF1) {
          isQRS = true;
          _SPKI = 0.125 * integVal + 0.875 * _SPKI;
          _SPKF = 0.125 * filtVal + 0.875 * _SPKF;
        }
        // Если не прошёл, но интервал превысил RR_MISSED_LIMIT – запускаем search‑back
        else if (_rrAvg1.isNotEmpty &&
            _rrAvg1.last > _RRMissedLimit) {
          final int searchBack = _performSearchBack(
            integrationSignal: integrationSignal,
            bandPassSignal: bandPassSignal,
            currentIntegIdx: integIdx,
            fs: fs,
          );
          if (searchBack != -1) {
            // Нашли пропущенный комплекс
            isQRS = true;
            // Обновляем SPK с коэффициентом 0.25 (для поиска назад)
            _SPKI = 0.25 * integrationSignal[searchBack] + 0.75 * _SPKI;
            _SPKF = 0.25 * bandPassSignal[searchBack] + 0.75 * _SPKF;
            // Добавляем найденный fiducial
            rPeaks.add(searchBack);
            _lastQRSIndex = searchBack;
            _updateThresholds();
            // Пропускаем текущий пик (он уже обработан)
            continue;
          }
        }

        // Если QRS подтверждён
        if (isQRS) {
          // ---- T-волновая дискриминация (если RR < 360 мс) ----
          if (_rrAvg1.isNotEmpty &&
              _rrAvg1.last < 0.36 &&
              _rrAvg1.last > 0.2) {
            final double currSlope = _computeMaxSlope(integrationSignal, integIdx, fs);
            if (currSlope < 0.5 * _prevSlope) {
              // Это T-волна – обновляем шумовой порог
              _NPKI = 0.125 * integVal + 0.875 * _NPKI;
              _NPKF = 0.125 * filtVal + 0.875 * _NPKF;
              _updateThresholds();
              continue; // не добавляем как QRS
            }
          }

          // Добавляем fiducial метку (пик интегрального сигнала)
          rPeaks.add(integIdx);
          _lastQRSIndex = integIdx;
          _prevSlope = _computeMaxSlope(integrationSignal, integIdx, fs);
          _updateThresholds();
        } else {
          // Пик не является QRS – обновляем шумовые оценки
          _NPKI = 0.125 * integVal + 0.875 * _NPKI;
          _NPKF = 0.125 * filtVal + 0.875 * _NPKF;
          _updateThresholds();
        }
      }

      return rPeaks;
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

    void _updateRRIntervals(List<_PeakPair> pairs, int currentIdx, int fs) {
      if (currentIdx < 1) return;
      final double rr = (pairs[currentIdx].integIdx - pairs[currentIdx - 1].integIdx) / fs;
      _rrAvg1.add(rr);
      if (_rrAvg1.length > 8) _rrAvg1.removeAt(0);

      // Обновление RR_AVERAGE2 (только RR в пределах 92–116% от текущего среднего)
      if (_rrAvg2.length >= 8) {
        final double avg2 = _rrAvg2.reduce((a, b) => a + b) / _rrAvg2.length;
        _RRLowLimit = 0.92 * avg2;
        _RRHighLimit = 1.16 * avg2;
        _RRMissedLimit = 1.66 * avg2;
      }
      if (_rrAvg2.isEmpty ||
          (rr >= _RRLowLimit && rr <= _RRHighLimit)) {
        _rrAvg2.add(rr);
        if (_rrAvg2.length > 8) _rrAvg2.removeAt(0);
      }
    }

    int _performSearchBack({
      required List<double> integrationSignal,
      required List<double> bandPassSignal,
      required int currentIntegIdx,
      required int fs,
    }) {
      // Ищем пик в интервале [currentIntegIdx - RR_MISSED_LIMIT*fs, currentIntegIdx]
      final int searchWindow = (_RRMissedLimit * fs).round();
      int start = (currentIntegIdx - searchWindow).clamp(0, integrationSignal.length - 1);
      int end = currentIntegIdx.clamp(0, integrationSignal.length - 1);

      int bestIntegIdx = -1;
      double maxVal = -double.infinity;
      // Ищем пик между нижним и верхним порогами (нижний порог используется)
      for (int i = start; i <= end; i++) {
        final double val = integrationSignal[i];
        if (val > _THRESHOLDI2 && val < _THRESHOLDI1 && val > maxVal) {
          maxVal = val;
          bestIntegIdx = i;
        }
      }
      if (bestIntegIdx == -1) return -1;

      // Проверяем соответствующий пик в фильтрованном сигнале (нижний порог)
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

      return bestIntegIdx; // возвращаем fiducial как пик интегрального сигнала
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

  /// Основной метод обработки сигнала (экземплярный)
  List<int> process(List<double> signal) {
    if (signal.isEmpty) return [];

    // Сброс всех внутренних состояний (буферы фильтров + пороги)
    reset();
    resetThresholds();

    // 1. Полосовая фильтрация (каскад ФНЧ + ФВЧ)
    final List<double> filtered = List.filled(signal.length, 0.0);
    for (int i = 0; i < signal.length; i++) {
      double low = _applyLowPassFilter(signal[i]);
      filtered[i] = _applyHighPassFilter(low) / 1152.0; // общий gain = 36*32
    }

    // 2. Производная
    final List<double> derivated = derivative(filtered, _sampleRate);

    // 3. Возведение в квадрат
    final List<double> squared = squaring(derivated);

    // 4. Скользящее интегрирование (окно 150 мс)
    final List<double> integrated = movingWindowIntegration(squared, _sampleRate);

    // 5. Детекция пиков
    final List<int> peaks = detectPeaks(
      ecgSignal: signal,          // не используется внутри, но оставлен для совместимости
      fs: _sampleRate.round(),
      integrationSignal: integrated,
      bandPassSignal: filtered,
    );

    // Возвращаем только индексы R-пиков
    return peaks;
  }
}

  // Вспомогательный класс для хранения пары индексов
  class _PeakPair {
    final int integIdx;
    final int filtIdx;
    _PeakPair({required this.integIdx, required this.filtIdx});
  }