import 'dart:io';
import 'dart:async';
import 'dart:convert';
import 'dart:math';

/// Глобальная переменная для пути к утилитам WFDB
/// По умолчанию пустая строка - утилиты ищутся в PATH
String wfdbBinPath = 'C:/Instruments/wfdb_utils/wfdb-software-package-10.6.2/build/bin';

class PanTompkinsQRS {
  List<double> bandPassFilter(List<double> signal) {
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
  }

  List<double> derivative(List<double> signal, double fs) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      result[index] = 0;

      if (index >= 1) {
        result[index] -= 2 * signal[index - 1];
      }

      if (index >= 2) {
        result[index] -= signal[index - 2];
      }

      if (index >= 2 && index <= signal.length - 2) {
        result[index] += 2 * signal[index + 1];
      }

      if (index >= 2 && index <= signal.length - 3) {
        result[index] += signal[index + 2];
      }

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
    List<double> bpass = bandPassFilter(inputSignal);
    List<double> der = derivative(bpass, fs.toDouble());
    List<double> sqr = squaring(der);
    List<double> mwin = movingWindowIntegration(sqr, fs.toDouble());
    List<int> peaks = detectPeaks(
        ecgSingal: signal,
        fs: fs,
        integration_signal: mwin,
        band_pass_signal: bpass);
    double heartRate = (60 * fs) / average(diff(peaks.sublist(1)));
    return (heartRate, peaks);
  }

  List<int> detectPeaks(
      {required List<double> ecgSingal,
      required int fs,
      required List<double> integration_signal,
      required List<double> band_pass_signal}) {
    List<int> possible_peaks = [];
    List<int> signal_peaks = [];
    List<int> r_peaks = [];
    double SPKI = 0;
    double SPKF = 0;
    double NPKI = 0;
    double NPKF = 0;
    List rr_avg_one = [];
    double THRESHOLDI1 = 0;
    double THRESHOLDF1 = 0;
    List<double> rr_avg_two = [];
    double THRESHOLDI2 = 0;
    double THRESHOLDF2 = 0;
    int is_T_found = 0;
    double current_slope = 0;
    double previous_slope = 0;
    int window = (0.15 * fs).round();
    List FM_peaks = [];
    List<double> integration_signal_smooth =
        convolution([...integration_signal]);
    List localDiff = diff([...integration_signal_smooth]);
    double RR_LOW_LIMIT = -999;
    double RR_HIGH_LIMIT = 999;
    double RR_MISSED_LIMIT = 0;

    for (int i = 1; i < localDiff.length; i++) {
      if (i - 1 > 2 * fs && localDiff[i - 1] > 0 && localDiff[i] < 0) {
        FM_peaks.add(i - 1);
      }
    }
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
      if (max_index != -1) {
        possible_peaks.add(max_index);
      }
      if (index == 0 || index > possible_peaks.length) {
        if (integration_signal[current_peak] >= THRESHOLDI1) {
          SPKI = 0.125 * integration_signal[current_peak] + 0.875 * SPKI;
          if (possible_peaks[index] > THRESHOLDF1) {
            SPKF = 0.125 * band_pass_signal[index] + 0.875 * SPKF;
            signal_peaks.add(possible_peaks[index]);
          } else {
            NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
          }
        } else if ((integration_signal[current_peak] > THRESHOLDI2 &&
                integration_signal[current_peak] < THRESHOLDI1) ||
            (integration_signal[current_peak] < THRESHOLDI2)) {
          NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
          NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
        }
      } else {
        List RRAVERAGE1 = divideList(
            diff(FM_peaks.sublist(max(0, index - 8), index + 1)), fs);
        double rr_one_mean = average(RRAVERAGE1);
        rr_avg_one.add(rr_one_mean);
        double limit_factor = rr_one_mean;

        if (index >= 8) {
          for (double RR in RRAVERAGE1) {
            if (RR > RR_LOW_LIMIT && RR < RR_HIGH_LIMIT) {
              rr_avg_two.add(RR);
              if (rr_avg_two.length == 9) {
                rr_avg_two.removeAt(0);
                limit_factor = average(rr_avg_two);
              }
            }
          }
          if (rr_avg_two.length == 8 || index < 8) {
            RR_LOW_LIMIT = 0.92 * limit_factor;
            RR_HIGH_LIMIT = 1.16 * limit_factor;
            RR_MISSED_LIMIT = 1.66 * limit_factor;
            RR_MISSED_LIMIT = 1.66 * limit_factor;
          }
          if (rr_avg_one[rr_avg_one.length - 1] < RR_LOW_LIMIT ||
              rr_avg_one[rr_avg_one.length - 1] > RR_MISSED_LIMIT) {
            THRESHOLDI1 = THRESHOLDI1 / 2;
            THRESHOLDF1 = THRESHOLDF1 / 2;
          }

          double curr_rr_interval = RRAVERAGE1[RRAVERAGE1.length - 1];
          int search_back_window = (curr_rr_interval * fs).round();
          if (curr_rr_interval > RR_MISSED_LIMIT) {
            left_limit = current_peak - search_back_window + 1;
            right_limit = current_peak + 1;
            int search_back_max_index = -1;
            max_value = -999999;
            for (int i = left_limit; i < right_limit; i++) {
              if (integration_signal[i] > THRESHOLDI1 &&
                  integration_signal[i] > max_value) {
                max_value = integration_signal[i];
                search_back_max_index = i;
              }
            }
            if (search_back_max_index != -1) {
              SPKI = 0.25 *
                      integration_signal[search_back_max_index < 0
                          ? integration_signal.length + search_back_max_index
                          : search_back_max_index] +
                  0.75 * SPKI;
              THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
              THRESHOLDI2 = 0.5 * THRESHOLDI1;

              left_limit = search_back_max_index - (0.15 * fs).round();
              right_limit = min(band_pass_signal.length, search_back_max_index);

              int search_back_max_index2 = -1;
              max_value = -999999;

              for (int i = left_limit; i < right_limit; i++) {
                if (band_pass_signal[i] > THRESHOLDF1 &&
                    band_pass_signal[i] > max_value) {
                  max_value = band_pass_signal[i];
                  search_back_max_index2 = i;
                }
              }

              if (search_back_max_index2 != -1 &&
                  band_pass_signal[search_back_max_index2 < 0
                      ? integration_signal.length + search_back_max_index2
                      : search_back_max_index2] > THRESHOLDF2) {
                SPKF = 0.25 *
                        band_pass_signal[search_back_max_index2 < 0
                            ? integration_signal.length + search_back_max_index2
                            : search_back_max_index2] +
                    0.75 * SPKF;
                THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
                THRESHOLDF2 = 0.5 * THRESHOLDF1;
                signal_peaks.add(search_back_max_index2);
              }
            }
          }
          if (integration_signal[current_peak] >= THRESHOLDI1) {
            if (curr_rr_interval > 0.20 &&
                curr_rr_interval < 0.36 &&
                index > 0) {
              current_slope = geMax(diff(integration_signal.sublist(
                  current_peak - (fs * 0.075).round(), current_peak + 1)));

              previous_slope = geMax(diff(integration_signal.sublist(
                  FM_peaks[index - 1] - (fs * 0.075).round(),
                  FM_peaks[index - 1] + 1)));
              if (current_slope < 0.5 * previous_slope) {
                NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
                is_T_found = 1;
              }
            }

            if (is_T_found == 0) {
              SPKI = 0.125 * integration_signal[current_peak] + 0.875 * SPKI;

              if (possible_peaks[index] > THRESHOLDF1) {
                SPKF = 0.125 * band_pass_signal[index] + 0.875 * SPKF;
                signal_peaks.add(possible_peaks[index]);
              } else {
                NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
              }
            }
          } else if ((integration_signal[current_peak] > THRESHOLDI1 &&
                  integration_signal[current_peak] < THRESHOLDI2) ||
              (integration_signal[current_peak] < THRESHOLDI1)) {
            NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
            NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
          }

          THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
          THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
          THRESHOLDI2 = 0.5 * THRESHOLDI1;
          THRESHOLDF2 = 0.5 * THRESHOLDF1;
          is_T_found = 0;
        }
      }
    }
    for (int i in unique(signal_peaks)) {
      int window = (0.2 * fs).round();
      int left_limit = i - window;
      double right_limit = min(i + window + 1, ecgSingal.length.toDouble());
      double max_value = -double.infinity;
      int max_index = -1;
      for (int i = left_limit; i < right_limit; i++) {
        if (ecgSingal[i] > max_value) {
          max_value = ecgSingal[i];
          max_index = i;
        }
      }
      r_peaks.add(max_index);
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