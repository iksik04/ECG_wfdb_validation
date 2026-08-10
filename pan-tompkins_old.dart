import 'dart:io';
import 'dart:async';
import 'dart:convert';
import 'dart:math';

/// Глобальная переменная для пути к утилитам WFDB
/// По умолчанию пустая строка - утилиты ищутся в PATH
String wfdbBinPath = 'C:/Instruments/wfdb_utils/wfdb-software-package-10.6.2/build/bin';

class PanTompkinsQRS {
  /*List<double> bandPassFilter(List<double> signal) {
    List<double> sig = List.from(signal);
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      sig[index] = signal[index];

      if (index >= 1) {
        sig[index] += 2 * sig[index - 1];
      }

      if (index >= 2) {
        sig[index] -= sig[index - 2];
      }

      if (index >= 6) {
        sig[index] -= 2 * signal[index - 6];
      }

      if (index >= 12) {
        sig[index] += signal[index - 12];
      }
    }

    result = List.from(sig);

    for (int index = 0; index < signal.length; index++) {
      result[index] = -1 * sig[index];

      if (index >= 1) {
        result[index] -= result[index - 1];
      }

      if (index >= 16) {
        result[index] += 32 * sig[index - 16];
      }

      if (index >= 32) {
        result[index] += sig[index - 32];
      }
    }

    double maxVal = max(result.reduce(max), -result.reduce(min));
    result = result.map((val) => val / maxVal).toList();

    return result;
  }*/

  List<double> derivative(List<double> signal, double fs) {
    List<double> result = List.filled(signal.length, 0.0);

    for (int index = 0; index < signal.length; index++) {
      result[index] = 0;

      if (index >= 2) result[index] -= signal[index - 2];
      if (index >= 1) result[index] -= 2 * signal[index - 1];
      if (index < signal.length - 2) result[index] += 2 * signal[index + 1];
      if (index < signal.length - 3) result[index] += signal[index + 2];

      result[index] = (result[index] * fs) / 8;
    }
    return result;
  }

  List<double> squaring(List<double> signal) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      result[index] = signal[index] * signal[index];
    }

    return result;
  }

  List<double> movingWindowIntegration(List<double> signal, double fs) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];
    int winSize = (0.150 * fs).round();
    double sum = 0;

    for (int j = 0; j < winSize; j++) {
      sum += signal[j] / winSize;
      result[j] = sum;
    }

    for (int index = winSize; index < signal.length; index++) {
      sum += signal[index] / winSize;
      sum -= signal[index - winSize] / winSize;
      result[index] = sum;
    }

    return result;
  }

  (double, List<int>) solve(List<double> signal, int fs) {
    List<double> inputSignal = List.from(signal);
    //List<double> bpass = bandPassFilter(inputSignal);
    List<double> der = derivative(inputSignal, fs.toDouble());
    List<double> sqr = squaring(der);
    List<double> mwin = movingWindowIntegration(sqr, fs.toDouble());
    List<int> peaks = detectPeaks(
        ecgSingal: signal,
        fs: fs,
        integration_signal: mwin,
        band_pass_signal: signal);
    double heartRate = (60 * fs) / average(diff(peaks.sublist(1)));
    return (heartRate, peaks);
  }

  // TODO: разобраться с фильтрами и тем, как они влияют на detectPeaks
  
  List<int> detectPeaks({
  required List<double> ecgSingal,
  required int fs,
  required List<double> integration_signal,
  required List<double> band_pass_signal,
}) {
  // ---- Объявление всех переменных (были удалены) ----
  double SPKI = 0.0;
  double NPKI = 0.0;
  double SPKF = 0.0;
  double NPKF = 0.0;
  double THRESHOLDI1 = 0.0;
  double THRESHOLDF1 = 0.0;
  double THRESHOLDI2 = 0.0;
  double THRESHOLDF2 = 0.0;
  
  List<int> possible_peaks = [];
  List<int> signal_peaks = [];
  List<int> r_peaks = [];

  // ---- Пункт 1: инициализация порогов по первым 2 секундам ----
  int initLen = (2 * fs).clamp(0, band_pass_signal.length);
  List<double> intInit = integration_signal.sublist(0, initLen);
  List<double> bpInit = band_pass_signal.sublist(0, initLen);

  List<int> localMaxIdx = [];
  for (int i = 1; i < initLen - 1; i++) {
    if (intInit[i] > intInit[i - 1] && intInit[i] > intInit[i + 1]) {
      localMaxIdx.add(i);
    }
  }
  if (localMaxIdx.isEmpty) {
    double meanInt = intInit.reduce((a, b) => a + b) / intInit.length;
    double meanBp = bpInit.reduce((a, b) => a + b) / bpInit.length;
    SPKI = meanInt * 0.5;
    NPKI = meanInt * 0.5;
    SPKF = meanBp * 0.5;
    NPKF = meanBp * 0.5;
  } else {
    List<double> intVals = localMaxIdx.map((i) => intInit[i]).toList();
    intVals.sort((a, b) => b.compareTo(a));
    List<double> bpVals = localMaxIdx.map((i) => bpInit[i]).toList();
    bpVals.sort((a, b) => b.compareTo(a));

    if (intVals.length >= 2) {
      SPKI = (intVals[0] + intVals[1]) / 2;
      NPKI = intVals.sublist(2).reduce((a, b) => a + b) / intVals.sublist(2).length;
      SPKF = (bpVals[0] + bpVals[1]) / 2;
      NPKF = bpVals.sublist(2).reduce((a, b) => a + b) / bpVals.sublist(2).length;
    } else {
      SPKI = intVals[0];
      NPKI = 0.0;
      SPKF = bpVals[0];
      NPKF = 0.0;
    }
  }

  THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
  THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
  THRESHOLDI2 = 0.5 * THRESHOLDI1;
  THRESHOLDF2 = 0.5 * THRESHOLDF1;

  int is_T_found = 0;
  double current_slope = 0;
  double previous_slope = 0;
  int window = (0.15 * fs).round();

  // ---- Пункт 4: исправленная свёртка (без сдвига) ----
  List<double> integration_signal_smooth = convolution([...integration_signal]);
  List localDiff = diff([...integration_signal_smooth]);

  List<int> FM_peaks = [];
  for (int i = 1; i < localDiff.length; i++) {
    if (i - 1 > 2 * fs && localDiff[i - 1] > 0 && localDiff[i] < 0) {
      FM_peaks.add(i - 1);
    }
  }

  // ---- Пункт 2: RR-буфер и лимиты ----
  List<double> rr_history = [];
  double RR_LOW_LIMIT = 0.2;
  double RR_HIGH_LIMIT = 0.8;
  double RR_MISSED_LIMIT = 1.2;

  // ---- Основной цикл по FM_peaks ----
  for (int index = 0; index < FM_peaks.length; index++) {
    int current_peak = FM_peaks[index];
    int left_limit = max(current_peak - window, 0);
    int right_limit = min(current_peak + window + 1, band_pass_signal.length);

    int max_index = -1;
    double max_value = -999999;
    for (int i = left_limit; i < right_limit; i++) {
      if (band_pass_signal[i] > max_value) {
        max_value = band_pass_signal[i];
        max_index = i;
      }
    }
    if (max_index == -1) continue;
    possible_peaks.add(max_index);

    double peak_int = integration_signal[current_peak];
    double peak_bp = band_pass_signal[max_index];

    // ---- Пункт 3: первый пик обрабатываем отдельно ----
    if (index == 0) {
      if (peak_int >= THRESHOLDI1) {
        SPKI = 0.125 * peak_int + 0.875 * SPKI;
        if (peak_bp > THRESHOLDF1) {
          SPKF = 0.125 * peak_bp + 0.875 * SPKF;
          signal_peaks.add(max_index);
        } else {
          NPKF = 0.125 * peak_bp + 0.875 * NPKF;
        }
      } else if (peak_int > THRESHOLDI2 && peak_int < THRESHOLDI1) {
        NPKI = 0.125 * peak_int + 0.875 * NPKI;
        NPKF = 0.125 * peak_bp + 0.875 * NPKF;
      }
      THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
      THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
      THRESHOLDI2 = 0.5 * THRESHOLDI1;
      THRESHOLDF2 = 0.5 * THRESHOLDF1;
      is_T_found = 0;
      continue;
    }

    // ---- Обработка остальных пиков ----
    double rr = (FM_peaks[index] - FM_peaks[index - 1]) / fs;
    rr_history.add(rr);
    if (rr_history.length > 8) rr_history.removeAt(0);

    if (rr_history.length >= 8) {
      double meanRR = rr_history.reduce((a, b) => a + b) / rr_history.length;
      RR_LOW_LIMIT = 0.92 * meanRR;
      RR_HIGH_LIMIT = 1.16 * meanRR;
      RR_MISSED_LIMIT = 1.66 * meanRR;
    }

    // ---- Search-back для пропущенных пиков ----
    if (rr > RR_MISSED_LIMIT) {
      int search_back_window = (rr * fs).round();
      left_limit = current_peak - search_back_window + 1;
      right_limit = current_peak + 1;
      int search_back_max_index = -1;
      double max_int = -999999;
      for (int i = left_limit; i < right_limit; i++) {
        if (i < 0 || i >= integration_signal.length) continue;
        if (integration_signal[i] > THRESHOLDI1 && integration_signal[i] > max_int) {
          max_int = integration_signal[i];
          search_back_max_index = i;
        }
      }
      if (search_back_max_index != -1) {
        SPKI = 0.25 * integration_signal[search_back_max_index] + 0.75 * SPKI;
        THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
        THRESHOLDI2 = 0.5 * THRESHOLDI1;

        int sb_left = search_back_max_index - (0.15 * fs).round();
        int sb_right = min(band_pass_signal.length, search_back_max_index);
        int sb_max_idx = -1;
        double max_bp = -999999;
        for (int i = sb_left; i < sb_right; i++) {
          if (i < 0 || i >= band_pass_signal.length) continue;
          if (band_pass_signal[i] > THRESHOLDF1 && band_pass_signal[i] > max_bp) {
            max_bp = band_pass_signal[i];
            sb_max_idx = i;
          }
        }
        if (sb_max_idx != -1 && band_pass_signal[sb_max_idx] > THRESHOLDF2) {
          SPKF = 0.25 * band_pass_signal[sb_max_idx] + 0.75 * SPKF;
          THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
          THRESHOLDF2 = 0.5 * THRESHOLDF1;
          signal_peaks.add(sb_max_idx);
        }
      }
    }

    // ---- Подавление T-зубцов ----
    if (rr > 0.20 && rr < 0.36 && index > 0) {
      int prev_peak = FM_peaks[index - 1];
      current_slope = geMax(diff(integration_signal.sublist(
          current_peak - (fs * 0.075).round(), current_peak + 1)));
      previous_slope = geMax(diff(integration_signal.sublist(
          prev_peak - (fs * 0.075).round(), prev_peak + 1)));
      if (current_slope < 0.5 * previous_slope) {
        NPKI = 0.125 * peak_int + 0.875 * NPKI;
        is_T_found = 1;
      }
    }

    // ---- Основное решение о пике ----
    if (is_T_found == 0) {
      if (peak_int >= THRESHOLDI1) {
        SPKI = 0.125 * peak_int + 0.875 * SPKI;
        if (peak_bp > THRESHOLDF1) {
          SPKF = 0.125 * peak_bp + 0.875 * SPKF;
          signal_peaks.add(max_index);
        } else {
          NPKF = 0.125 * peak_bp + 0.875 * NPKF;
        }
      } else if (peak_int > THRESHOLDI2 && peak_int < THRESHOLDI1) {
        NPKI = 0.125 * peak_int + 0.875 * NPKI;
        NPKF = 0.125 * peak_bp + 0.875 * NPKF;
      }
    }

    // ---- Обновление порогов ----
    THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
    THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
    THRESHOLDI2 = 0.5 * THRESHOLDI1;
    THRESHOLDF2 = 0.5 * THRESHOLDF1;
    is_T_found = 0;
  }

  // ---- Формирование финальных R-пиков ----
  for (int i in unique(signal_peaks)) {
    int w = (0.2 * fs).round();
    int left_limit = i - w;
    int right_limit = min(i + w + 1, ecgSingal.length);
    double max_val = -double.infinity;
    int max_idx = -1;
    for (int j = left_limit; j < right_limit; j++) {
      if (j < 0) continue;
      if (ecgSingal[j] > max_val) {
        max_val = ecgSingal[j];
        max_idx = j;
      }
    }
    if (max_idx != -1) r_peaks.add(max_idx);
  }

  return r_peaks;
}
}
/// Вспомогательные функции
List<double> convolution(List<double> signal) {
  List<double> sig1 = [...signal];
  List<double> sig2 = [for (int i = 0; i < 20; i++) 0.05];
  List<double> conv = [for (int i = 0; i < (sig1.length - sig2.length); i++) 0];
  for (int l = 0; l < conv.length; l++) {
    for (int i = 0; i < sig2.length; i++) {
      conv[l] += sig1[l - i + sig2.length] * sig2[i];
    }
  }

  List<double> d = [for (int i = 0; i < signal.length; i++) 0];
  for (int l = 0; l < conv.length; l++) {
    d[l + 11] = conv[l];
  }
  return d;
}

double geMax(List signal) {
  double largestGeekValue = signal[0];

  for (var i = 0; i < signal.length; i++) {
    if (signal[i] > largestGeekValue) {
      largestGeekValue = signal[i];
    }
  }
  return largestGeekValue;
}

List diff(List signal) {
  List out = [];
  if (signal.length < 2) {
    return [];
  }
  for (int i = 0; i < signal.length - 1; i++) {
    out.add(signal[i + 1] - signal[i]);
  }
  return out;
}

List unique(List arr) {
  arr.sort();
  List unique_list = [];

  var last_added;

  for (var element in arr) {
    if (element != last_added) {
      unique_list.add(element);
      last_added = element;
    }
  }
  return unique_list;
}

List divideList(List signal, int divider) {
  List out = [];
  for (int i = 0; i < signal.length; i++) {
    out.add(signal[i] / divider);
  }
  return out;
}

double average(List signal) {
  double summ = 0;
  for (int i = 0; i < signal.length; i++) {
    summ += signal[i];
  }
  return summ / signal.length;
}

/// Получить полный путь к утилите WFDB
String getWfdbCommand(String command) {
  if (wfdbBinPath.isEmpty) {
    return command;
  }
  // Добавляем разделитель пути, если его нет в конце
  String path = wfdbBinPath;
  if (!path.endsWith(Platform.pathSeparator)) {
    path += Platform.pathSeparator;
  }
  // На Windows добавляем .exe если нужно
  if (Platform.isWindows && !command.endsWith('.exe')) {
    return path + command + '.exe';
  }
  return path + command;
}

/// Чтение частоты дискретизации из файла .hea
Future<double> getSampleRate(String filePath) async {
  try {
    // Нормализуем путь для Windows
    String normalizedPath = filePath.replaceAll('/', '\\');
    if (normalizedPath.endsWith('.hea')) {
      normalizedPath = normalizedPath.substring(0, normalizedPath.length - 4);
    }
    
    final heaFile = File('$normalizedPath.hea');
    if (!await heaFile.exists()) {
      print('Предупреждение: файл .hea не найден, используется частота по умолчанию 360');
      return 360.0;
    }

    final content = await heaFile.readAsString();
    final lines = content.split('\n');
    
    for (final line in lines) {
      if (line.trim().isEmpty) continue;
      
      final parts = line.trim().split(RegExp(r'\s+'));
      if (parts.length >= 3) {
        final sampleRate = double.tryParse(parts[2]);
        if (sampleRate != null && sampleRate > 0) {
          return sampleRate;
        }
      }
      break;
    }

    print('Предупреждение: не удалось определить частоту дискретизации, используется 360');
    return 360.0;
  } catch (e) {
    print('Ошибка чтения файла .hea: $e');
    return 360.0;
  }
}

/// Загрузка данных ЭКГ с помощью rdsamp
Future<List<double>> loadECGDataWithRDSamp(String folderPath, String recordName, int channel) async {
  try {
    // Убеждаемся, что файл физически существует по указанному пути
    final datFile = File('$folderPath\\$recordName.dat');
    if (!await datFile.exists()) {
      print('Ошибка: файл .dat не найден: ${datFile.path}');
      return [];
    }

    final rdsampCmd = getWfdbCommand('rdsamp');
    print('Запуск rdsamp для записи: $recordName в папке: $folderPath');

    // Используем Process.run с установкой окружения и рабочей директории
    final result = await Process.run(
      rdsampCmd,
      // ВАЖНО: передаем только имя записи, а не полный путь
      ['-r', recordName, '-f', '0', '-t', 'end', '-p', '-v'],
      runInShell: true,
      workingDirectory: folderPath, // Выполняем команду прямо в папке с БД
      environment: {'WFDB': '.'},    // Принудительно заставляем искать файлы локально
    );

    if (result.exitCode != 0) {
      print('Ошибка выполнения rdsamp (код ${result.exitCode}): ${result.stderr}');
      return [];
    }

    return _parseRDSampOutput(result.stdout.toString(), channel);
  } catch (e) {
    print('Ошибка выполнения rdsamp: $e');
    return [];
  }
}

List<double> _parseRDSampOutput(String output, int channel) {
  final lines = output.split('\n')
      .where((line) => line.trim().isNotEmpty)
      .toList();

  if (lines.isEmpty) {
    print('Вывод rdsamp пуст');
    return [];
  }

  final signal = <double>[];

  for (final line in lines) {
    try {
      final parts = line.trim().split(RegExp(r'\s+'));
      if (parts.length < channel + 2) continue;

      final value = double.parse(parts[channel + 1]);

      if (value.isFinite) {
        signal.add(value);
      }
    } catch (e) {
      continue;
    }
  }
  return signal;
}

/// Запись пиков в аннотационный файл с помощью wrann
Future<void> writePeaksWithWRAnn(String folderPath, String recordName, List<int> peaks, int fs) async {
  try {
    final datFile = File('$folderPath\\$recordName.dat');
    if (!await datFile.exists()) {
      print('Ошибка: файл .dat не найден: ${datFile.path}');
      return;
    }

    final content = StringBuffer();
    for (int peak in peaks) {
      double timeInSeconds = peak / fs;
      int minutes = timeInSeconds ~/ 60;
      double seconds = timeInSeconds % 60;
      int wholeSeconds = seconds.floor();
      int milliseconds = ((seconds - wholeSeconds) * 1000).round();
      
      String formattedTime = '${minutes.toString().padLeft(1, '0')}:${wholeSeconds.toString().padLeft(2, '0')}.${milliseconds.toString().padLeft(3, '0')}';
      content.writeln('$formattedTime       $peak     N');
    }
    
    final wrannCmd = getWfdbCommand('wrann');
    
    // Используем Process.start с установкой окружения и рабочей директории
    final process = await Process.start(
      wrannCmd,
      // ВАЖНО: передаем только имя записи
      ['-r', recordName, '-a', 'gqrs'],
      mode: ProcessStartMode.normal,
      workingDirectory: folderPath, // Переходим в папку с БД
      environment: {'WFDB': '.'},    // Задаем переменную окружения
    );

    process.stdin.write(content.toString());
    await process.stdin.close();

    final output = await process.stdout.transform(utf8.decoder).join();
    final stderr = await process.stderr.transform(utf8.decoder).join();

    final exitCode = await process.exitCode;

    if (exitCode != 0) {
      print('Ошибка wrann (код $exitCode): $stderr');
      if (stderr.isNotEmpty) print('stderr: $stderr');
      if (output.isNotEmpty) print('stdout: $output');
      return;
    }
    print('Аннотации gqrs успешно записаны для $recordName');
  } catch (e) {
    print('Ошибка записи аннотационного файла: $e');
  }
}

/// Проверка доступности rdsamp
Future<bool> isRDSampAvailable() async {
  try {
    return true;
  } catch (e) {
    try {
      if (Platform.isWindows) {
        final result = await Process.run('where', ['rdsamp']);
        return result.exitCode == 0;
      } else {
        final result = await Process.run('which', ['rdsamp']);
        return result.exitCode == 0;
      }
    } catch (e2) {
      return false;
    }
  }
}

/// Проверка доступности wrann
Future<bool> isWRAnnAvailable() async {
  try {
    return true;
  } catch (e) {
    try {
      if (Platform.isWindows) {
        final result = await Process.run('where', ['wrann']);
        return result.exitCode == 0;
      } else {
        final result = await Process.run('which', ['wrann']);
        return result.exitCode == 0;
      }
    } catch (e2) {
      return false;
    }
  }
}

/// Обработка одной записи
Future<void> processRecording(String folderPath, String recordNumber, int channel) async {
  // Нормализуем пути
  String normalizedFolder = folderPath.replaceAll('/', '\\');
  // Убираем возможный trailing slash
  if (normalizedFolder.endsWith('\\')) {
    normalizedFolder = normalizedFolder.substring(0, normalizedFolder.length - 1);
  }
  
  print('Обработка: $normalizedFolder, запись $recordNumber, канал $channel');
  
  // Проверяем существование .dat файла
  final datFile = File('$normalizedFolder\\$recordNumber.dat');
  if (!await datFile.exists()) {
    print('Ошибка: файл $normalizedFolder\\$recordNumber.dat не найден');
    return;
  }

  // Получаем частоту дискретизации
  final sampleRate = await getSampleRate('$normalizedFolder\\$recordNumber');
  final fs = sampleRate.round();

  // Загружаем данные через rdsamp
  final ecgData = await loadECGDataWithRDSamp(folderPath, recordNumber, channel);
  if (ecgData.isEmpty) {
    print('Ошибка: данные не загружены');
    return;
  }

  print('Загружено ${ecgData.length} отсчётов, частота $fs Гц');

  // Расчитываем коэффициенты фильтров
  calculateFilterCoefficients(sampleRate);

  // Применяем фильтры
  List<double> filteredData = [for (double i in ecgData) applyHighPassFilter(i)];
  filteredData = [for (double i in filteredData) applyLowPassFilter(i)];

  // Детектируем пики
  final detector = PanTompkinsQRS();
  late List<int> peaks;
  late double heartRate;
  (heartRate, peaks) = detector.solve(filteredData, fs);
  
  print('Результаты: ${heartRate.toStringAsFixed(2)} BPM, ${peaks.length} пиков');

  // Сохраняем пики в .gqrs аннотацию
  await writePeaksWithWRAnn(folderPath, recordNumber, peaks, fs);

  // Сбрасываем состояние фильтров
  hprevFilterd = 0.0;
  hprevUnFiltered = 0.0;
  hprevprevUnfiltered = 0.0;
  hprevprevFilterd = 0.0;
  lprevFilterd = 0.0;
  lprevUnFiltered = 0.0;
  lprevprevUnfiltered = 0.0;
  lprevprevFilterd = 0.0;

  bhp = [];
  ahp = [];
  blp = [];
  alp = [];
}

void main(List<String> args) async {
  // Проверяем аргументы командной строки
  if (args.length < 3) {
    print('Использование: dart pan-tompkins.dart <путь_к_папке> <номер_записи> <канал>');
    print('Пример: dart pan-tompkins.dart ./assets/ECG_DB/AHADB 1201 1');
    print('');
    print('Аргументы:');
    print('  путь_к_папке   - Путь к папке с файлами записи');
    print('  номер_записи   - Номер записи (например, 1201)');
    print('  канал          - Номер канала (0 или 1)');
    print('');
    print('Программа:');
    print('  - Загружает данные через rdsamp');
    print('  - Детектирует R-пики алгоритмом Pan-Tompkins');
    print('  - Сохраняет пики в .gqrs аннотацию с помощью wrann');
    print('');
    print('Для указания пути к утилитам WFDB установите переменную wfdbBinPath');
    print('Например: wfdbBinPath = "C:/WFDB/bin/";');
    exit(1);
  }

  final folderPath = args[0];
  final recordNumber = args[1];
  final channel = int.parse(args[2]);

  await processRecording(folderPath, recordNumber, channel);
}

/// Переменные для фильтра высоких частот
double hprevFilterd = 0.0;
double hprevUnFiltered = 0.0;
double hprevprevUnfiltered = 0.0;
double hprevprevFilterd = 0.0;

List<double> bhp = [];
List<double> ahp = [];

/// Переменные для фильтра низких частот
double lprevFilterd = 0.0;
double lprevUnFiltered = 0.0;
double lprevprevUnfiltered = 0.0;
double lprevprevFilterd = 0.0;

List<double> blp = [];
List<double> alp = [];

void calculateFilterCoefficients(double sampleRate) {
  // Округляем до ближайшего целого, чтобы избежать ошибок округления double
  final int rate = sampleRate.round();

  if (rate == 125) {
    // --- ФНЧ (Low-Pass) 125 Гц ---
    blp = [0.35034638, 0.70069276, 0.35034638];
    alp = [-0.22115344, -0.18023207];

    // --- ФВЧ (High-Pass) 125 Гц ---
    bhp = [0.96851735, -1.93703469, 0.96851735];
    ahp = [1.93604329, -0.9380261];
    
  } else if (rate == 250) {
    // --- ФНЧ (Low-Pass) 250 Гц ---
    blp = [0.11735104, 0.23470207, 0.11735104];
    alp = [0.82523238, -0.29463653];

    // --- ФВЧ (High-Pass) 250 Гц ---
    bhp = [0.98413284, -1.96826569, 0.98413284];
    ahp = [1.96801391, -0.96851747];
    
  } else if (rate == 360) {
    // --- ФНЧ (Low-Pass) 360 Гц ---
    blp = [0.06433216, 0.12866431, 0.06433216];
    alp = [1.16557175, -0.42290037];

    // --- ФВЧ (High-Pass) 360 Гц ---
    bhp = [0.98895425, -1.9779085, 0.98895425];
    ahp = [1.97778648, -0.97803051];
  }
}

/// Фильтр низких частот
double applyLowPassFilter(double val) {
  double y = blp[0] * val +
      alp[0] * lprevFilterd +
      blp[1] * lprevUnFiltered +
      alp[1] * lprevprevFilterd +
      blp[2] * lprevprevUnfiltered;
  lprevprevFilterd = lprevFilterd;
  lprevFilterd = y;
  lprevprevUnfiltered = lprevUnFiltered;
  lprevUnFiltered = val;
  return y;
}

/// Фильтр высоких частот
double applyHighPassFilter(double val) {
  double y = bhp[0] * val +
      ahp[0] * hprevFilterd +
      bhp[1] * hprevUnFiltered +
      ahp[1] * hprevprevFilterd +
      bhp[2] * hprevprevUnfiltered;
  hprevprevFilterd = hprevFilterd;
  hprevFilterd = y;
  hprevprevUnfiltered = hprevUnFiltered;
  hprevUnFiltered = val;
  return y;
}