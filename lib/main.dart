// bin/main.dart
import 'dart:io';
import '../lib/wfdb_utils.dart';

// ============================================================
// Главная функция (точка входа)
// ============================================================
void main(List<String> args) async {
  // Проверяем аргументы командной строки
  if (args.length < 3) {
    print('Использование: dart main.dart <путь_к_папке> <номер_записи> <канал>');
    print('Пример: dart main.dart ./assets/ECG_DB/AHADB 1201 1');
    print('');
    print('Аргументы:');
    print('  путь_к_папке   - Путь к папке с файлами записи');
    print('  номер_записи   - Номер записи (например, 1201)');
    print('  канал          - Номер канала (0 или 1)');
    print('');
    print('Программа:');
    print('  - Загружает данные через rdsamp');
    print('  - Детектирует R-пики алгоритмом Pan-Tompkins (обновлённая версия)');
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