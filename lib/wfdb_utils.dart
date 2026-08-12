// lib/wfdb_utils.dart
import 'dart:io';
import 'dart:convert';
import 'dart:async';
import 'dart:math';

import '../lib/wqrs.dart';

// ============================================================
// Глобальная переменная для пути к утилитам WFDB
// ============================================================
String wfdbBinPath = 'C:/Instruments/wfdb_utils/wfdb-software-package-10.6.2/build/bin';

/// Получить полный путь к утилите WFDB
String getWfdbCommand(String command) {
  if (wfdbBinPath.isEmpty) {
    return command;
  }
  String path = wfdbBinPath;
  if (!path.endsWith(Platform.pathSeparator)) {
    path += Platform.pathSeparator;
  }
  if (Platform.isWindows && !command.endsWith('.exe')) {
    return path + command + '.exe';
  }
  return path + command;
}

/// Чтение коэффициента усиления из файла .hea
Future<double> getGain(String filePath) async {
  try {
    String normalizedPath = filePath.replaceAll('/', '\\');
    if (normalizedPath.endsWith('.hea')) {
      normalizedPath = normalizedPath.substring(0, normalizedPath.length - 4);
    }

    final heaFile = File('$normalizedPath.hea');
    if (!await heaFile.exists()) {
      print('Предупреждение: файл .hea не найден, используется усиление по умолчанию 200');
      return 200.0;
    }

    final content = await heaFile.readAsString();
    final lines = content.split('\n');

    for (int i = 0; i < lines.length; i++) {
      final line = lines[i].trim();
      if (line.isEmpty) continue;

      final parts = line.split(RegExp(r'\s+'));
      // Формат: имя_сигнала частота количество_бит ADC_разрешение [физический_минимум физический_максимум] [единицы]
      // Пример: 1201 2 250 16 0.0 100.0 mV
      if (parts.length >= 5) {
        final adcResolution = int.tryParse(parts[3]);
        final physicalMin = double.tryParse(parts[4]);
        final physicalMax = double.tryParse(parts[5]);
        
        if (adcResolution != null && physicalMin != null && physicalMax != null && physicalMax > physicalMin) {
          // Gain = (максимальное АЦП) / (физический диапазон)
          // Для 16-битного АЦП: max ADC = 2^16 = 65536
          // Но часто используется ±32768 для 16-битных данных
          double maxAdc = pow(2, adcResolution).toDouble() / 2;
          double physicalRange = physicalMax - physicalMin;
          double gain = maxAdc / physicalRange;
          print('Найден gain из .hea: $gain ADC/ед. (resolution=${adcResolution}bit, range=${physicalRange}ед.)');
          return gain;
        }
      }
    }

    print('Предупреждение: не удалось определить усиление, используется 200');
    return 200.0;
  } catch (e) {
    print('Ошибка чтения файла .hea: $e');
    return 200.0;
  }
}

/// Чтение частоты дискретизации из файла .hea
Future<double> getSampleRate(String filePath) async {
  try {
    String normalizedPath = filePath.replaceAll('/', '\\');
    if (normalizedPath.endsWith('.hea')) {
      normalizedPath = normalizedPath.substring(0, normalizedPath.length - 4);
    }

    final heaFile = File('$normalizedPath.hea');
    if (!await heaFile.exists()) {
      print('Предупреждение: файл .hea не найден, используется частота по умолчанию 250');
      return 250.0;
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

    print('Предупреждение: не удалось определить частоту дискретизации, используется 250');
    return 250.0;
  } catch (e) {
    print('Ошибка чтения файла .hea: $e');
    return 250.0;
  }
}

/// Загрузка сырых данных ЭКГ с помощью rdsamp (без флага -p)
Future<List<int>> loadRawECGDataWithRDSamp(String folderPath, String recordName, int channel) async {
  try {
    final datFile = File('$folderPath\\$recordName.dat');
    if (!await datFile.exists()) {
      print('Ошибка: файл .dat не найден: ${datFile.path}');
      return [];
    }

    final rdsampCmd = getWfdbCommand('rdsamp');
    print('Запуск rdsamp для записи: $recordName в папке: $folderPath');

    // БЕЗ флага -p - возвращаем сырые АЦП-значения
    final result = await Process.run(
      rdsampCmd,
      ['-r', recordName, '-f', '0', '-t', 'end'],  // Убрали -p
      runInShell: true,
      workingDirectory: folderPath,
      environment: {'WFDB': '.'},
    );

    if (result.exitCode != 0) {
      print('Ошибка выполнения rdsamp (код ${result.exitCode}): ${result.stderr}');
      return [];
    }

    return _parseRawRDSampOutput(result.stdout.toString(), channel);
  } catch (e) {
    print('Ошибка выполнения rdsamp: $e');
    return [];
  }
}

/// Парсинг сырых АЦП-данных из rdsamp
List<int> _parseRawRDSampOutput(String output, int channel) {
  final lines = output.split('\n')
      .where((line) => line.trim().isNotEmpty)
      .toList();

  if (lines.isEmpty) {
    print('Вывод rdsamp пуст');
    return [];
  }

  final signal = <int>[];

  for (final line in lines) {
    try {
      final parts = line.trim().split(RegExp(r'\s+'));
      if (parts.length < channel + 2) continue;

      final value = int.parse(parts[channel + 1]);
      signal.add(value);
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

      String formattedTime =
          '${minutes.toString().padLeft(1, '0')}:${wholeSeconds.toString().padLeft(2, '0')}.${milliseconds.toString().padLeft(3, '0')}';
      content.writeln('$formattedTime       $peak     N');
    }

    final wrannCmd = getWfdbCommand('wrann');

    final process = await Process.start(
      wrannCmd,
      ['-r', recordName, '-a', 'gqrs'],
      mode: ProcessStartMode.normal,
      workingDirectory: folderPath,
      environment: {'WFDB': '.'},
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

/// Обработка одной записи
Future<void> processRecording(String folderPath, String recordNumber, int channel) async {
  String normalizedFolder = folderPath.replaceAll('/', '\\');
  if (normalizedFolder.endsWith('\\')) {
    normalizedFolder = normalizedFolder.substring(0, normalizedFolder.length - 1);
  }

  print('Обработка: $normalizedFolder, запись $recordNumber, канал $channel');

  final datFile = File('$normalizedFolder\\$recordNumber.dat');
  if (!await datFile.exists()) {
    print('Ошибка: файл $normalizedFolder\\$recordNumber.dat не найден');
    return;
  }

  // Получаем параметры из .hea файла
  final sampleRate = await getSampleRate('$normalizedFolder\\$recordNumber');
  final gain = await getGain('$normalizedFolder\\$recordNumber');
  final fs = sampleRate.round();

  // Загружаем сырые АЦП-данные через rdsamp (без -p)
  final rawData = await loadRawECGDataWithRDSamp(folderPath, recordNumber, channel);
  if (rawData.isEmpty) {
    print('Ошибка: данные не загружены');
    return;
  }

  print('Загружено ${rawData.length} отсчётов, частота $fs Гц, усиление ${gain.toStringAsFixed(1)} ADC/ед.');

  // Создаём детектор для работы с сырыми АЦП-данными
  final detector = WQRSDetector(
    sampleRate: sampleRate,
    gain: gain,  // Используем реальное усиление из заголовка
    powerLineFreq: 50,
    minThreshold: 100,
    eyeClosingPeriod: 0.25,
    maxQRSWidth: 0.13,
    ndp: 2.5,
    debugMode: true,
    inputIsRaw: true,  // Указываем, что входные данные - сырые АЦП
  );

  List<int> peaks = detector.detect(rawData);

  print('Результат: ${peaks.length} пиков');

  await writePeaksWithWRAnn(folderPath, recordNumber, peaks, fs);
  
  var debugInfo = detector.getDebugInfo();
  print('\n=== DEBUG INFO ===');
  debugInfo.forEach((key, value) {
    if (value is double) {
      print('$key: ${value.toStringAsFixed(2)}');
    } else {
      print('$key: $value');
    }
  });
}