// lib/wfdb_utils.dart
import 'dart:io';
import 'dart:convert';
import 'dart:async';

import '../lib/pan-tompkins.dart';

// ============================================================
// Глобальная переменная для пути к утилитам WFDB
// ============================================================
String wfdbBinPath = 'C:/Instruments/wfdb_utils/wfdb-software-package-10.6.2/build/bin';

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
          print(sampleRate);
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
      ['-r', recordName, '-f', '0', '-t', 'end', '-p', '-v'],
      runInShell: true,
      workingDirectory: folderPath,
      environment: {'WFDB': '.'},
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
  // Нормализуем пути
  String normalizedFolder = folderPath.replaceAll('/', '\\');
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

  // Создаём детектор с нужной частотой дискретизации
  final detector = PanTompkinsQRS(sampleRate: sampleRate);

  // Детектируем пики
  final List<int> peaks = detector.detectRPeaks(ecgData);

  print('Результат: ${peaks.length} пиков');

  // Сохраняем пики в .gqrs аннотацию
  await writePeaksWithWRAnn(folderPath, recordNumber, peaks, fs);
}