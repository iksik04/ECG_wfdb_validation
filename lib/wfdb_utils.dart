// lib/wfdb_utils.dart
import 'dart:io';
import 'dart:convert';
import 'dart:async';

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

/// Информация о сигнале из файла .hea
class SignalInfo {
  final double sampleRate;
  final double gain;
  final int bits;
  final int adcZero;
  final int baseline;
  final String units;
  final String signalName;

  SignalInfo({
    required this.sampleRate,
    required this.gain,
    required this.bits,
    required this.adcZero,
    required this.baseline,
    required this.units,
    required this.signalName,
  });

  @override
  String toString() {
    return 'SignalInfo(sampleRate: $sampleRate, gain: $gain, bits: $bits, '
        'adcZero: $adcZero, baseline: $baseline, units: $units, signalName: $signalName)';
  }
}

/// Чтение информации о сигнале из файла .hea
Future<SignalInfo?> getSignalInfo(String filePath, int channel) async {
  try {
    // Нормализуем путь для Windows
    String normalizedPath = filePath.replaceAll('/', '\\');
    if (normalizedPath.endsWith('.hea')) {
      normalizedPath = normalizedPath.substring(0, normalizedPath.length - 4);
    }

    final heaFile = File('$normalizedPath.hea');
    if (!await heaFile.exists()) {
      print('Предупреждение: файл .hea не найден');
      return null;
    }

    final content = await heaFile.readAsString();
    final lines = content.split('\n');

    if (lines.isEmpty) {
      print('Предупреждение: файл .hea пуст');
      return null;
    }

    // Парсим первую строку - общая информация
    final firstParts = lines[0].trim().split(RegExp(r'\s+'));
    if (firstParts.length < 3) {
      print('Предупреждение: неверный формат .hea файла');
      return null;
    }

    final sampleRate = double.tryParse(firstParts[2]);
    if (sampleRate == null || sampleRate <= 0) {
      print('Предупреждение: неверная частота дискретизации');
      return null;
    }

    // Ищем строку для указанного канала (строка 1 + channel)
    if (lines.length <= channel + 1) {
      print('Предупреждение: канал $channel не найден в .hea файле');
      return null;
    }

    final channelLine = lines[channel + 1].trim();
    final channelParts = channelLine.split(RegExp(r'\s+'));

    if (channelParts.length < 6) {
      print('Предупреждение: неверный формат строки канала');
      return null;
    }

    // Формат строки канала:
    // fileName format gain bits adc_zero baseline [units] [signalName]
    // где:
    // fileName - имя файла данных (может быть пустым)
    // format - формат данных (16, 32, 64, или другие)
    // gain - усиление в мкВ/единицу АЦП
    // bits - количество бит
    // adc_zero - нулевое значение АЦП
    // baseline - базовое значение (обычно adc_zero)
    // units - единицы измерения (опционально)
    // signalName - имя сигнала (опционально)

    int index = 0;
    // Пропускаем имя файла если оно есть
    if (channelParts[0].isNotEmpty && !channelParts[0].contains(RegExp(r'^[0-9]+$'))) {
      index = 1;
    }

    // format
    final format = int.tryParse(channelParts[index]);
    if (format == null) {
      print('Предупреждение: неверный формат данных');
      return null;
    }
    index++;

    // gain
    final gain = double.tryParse(channelParts[index]);
    if (gain == null || gain == 0) {
      print('Предупреждение: неверное усиление');
      return null;
    }
    index++;

    // bits
    final bits = int.tryParse(channelParts[index]);
    if (bits == null) {
      print('Предупреждение: неверное количество бит');
      return null;
    }
    index++;

    // adc_zero
    final adcZero = int.tryParse(channelParts[index]);
    if (adcZero == null) {
      print('Предупреждение: неверное нулевое значение АЦП');
      return null;
    }
    index++;

    // baseline
    final baseline = int.tryParse(channelParts[index]);
    if (baseline == null) {
      print('Предупреждение: неверное базовое значение');
      return null;
    }
    index++;

    // units (опционально)
    String units = 'mV';
    if (channelParts.length > index) {
      units = channelParts[index];
      index++;
    }

    // signalName (опционально)
    String signalName = '';
    if (channelParts.length > index) {
      signalName = channelParts.sublist(index).join(' ');
    }

    return SignalInfo(
      sampleRate: sampleRate,
      gain: gain,
      bits: bits,
      adcZero: adcZero,
      baseline: baseline,
      units: units,
      signalName: signalName,
    );
  } catch (e) {
    print('Ошибка чтения файла .hea: $e');
    return null;
  }
}

/// Чтение частоты дискретизации из файла .hea (для обратной совместимости)
Future<double> getSampleRate(String filePath) async {
  final info = await getSignalInfo(filePath, 0);
  return info?.sampleRate ?? 360.0;
}

/// Загрузка данных ЭКГ с помощью rdsamp (без флага -p для сырых значений АЦП)
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
    // БЕЗ флага -p для получения сырых значений АЦП
    final result = await Process.run(
      rdsampCmd,
      ['-r', recordName, '-f', '0', '-t', 'end', '-v'],
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

      // При использовании rdsamp без -p, значения - это сырые значения АЦП
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

/// Загрузка данных ЭКГ с преобразованием в физические единицы
/// Использует gain из .hea для преобразования АЦП -> физические единицы
Future<List<double>> loadECGDataPhysical(String folderPath, String recordName, int channel) async {
  try {
    // Получаем информацию о сигнале
    final signalInfo = await getSignalInfo('$folderPath\\$recordName', channel);
    if (signalInfo == null) {
      print('Ошибка: не удалось получить информацию о сигнале');
      return [];
    }

    // Загружаем сырые данные
    final rawData = await loadECGDataWithRDSamp(folderPath, recordName, channel);
    if (rawData.isEmpty) {
      return [];
    }

    // Преобразуем в физические единицы (например, мкВ)
    // Формула: physical = (raw - baseline) / gain
    final physicalData = rawData.map((raw) {
      return (raw - signalInfo.baseline) / signalInfo.gain;
    }).toList();

    return physicalData;
  } catch (e) {
    print('Ошибка загрузки физических данных: $e');
    return [];
  }
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
      ['-r', recordName, '-a', 'wqrs'],
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
    print('Аннотации wqrs успешно записаны для $recordName');
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

  // Получаем информацию о сигнале
  final signalInfo = await getSignalInfo('$normalizedFolder\\$recordNumber', channel);
  if (signalInfo == null) {
    print('Ошибка: не удалось получить информацию о сигнале');
    return;
  }

  final fs = signalInfo.sampleRate.round();
  print('Частота дискретизации: $fs Гц, усиление: ${signalInfo.gain} мкВ/единицу АЦП');

  // Загружаем сырые данные (АЦП) через rdsamp без флага -p
  final ecgData = await loadECGDataWithRDSamp(folderPath, recordNumber, channel);
  if (ecgData.isEmpty) {
    print('Ошибка: данные не загружены');
    return;
  }

  print('Загружено ${ecgData.length} отсчётов (сырые значения АЦП)');

  // Создаём детектор с нужной частотой дискретизации
  // Передаем усиление в детектор для корректного расчета lfsc
  final detector = WQRS(
    sampleRate: signalInfo.sampleRate,
    signal: channel,
    gain: signalInfo.gain, // Теперь передается double
  );

  // Детектируем QRS комплексы
  final detections = detector.process(ecgData);
  
  // Извлекаем позиции QRS onset
  final List<int> peaks = detections.map((d) => d['time'] as int).toList();

  print('Найдено QRS комплексов: ${peaks.length}');
  
  // Сохраняем аннотации
  await writePeaksWithWRAnn(folderPath, recordNumber, peaks, fs);
}