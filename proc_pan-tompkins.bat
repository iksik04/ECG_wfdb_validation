@echo off
setlocal

set SCRIPT_NAME=pan-tompkins.dart
set DB_PATH=C:\Users\Public\Work\ECG_DB\AHADB

echo Searching for records in folder %DB_PATH%...
echo.

set WFDB=.

set COUNT=0

for %%f in ("%DB_PATH%\*.dat") do (
    echo Found record: %%~nf
    echo Running: dart %SCRIPT_NAME% %DB_PATH% %%~nf 0
    dart %SCRIPT_NAME% %DB_PATH% %%~nf 0
    echo.
    set /a COUNT+=1
)

echo.
echo Done! Processed records: %COUNT%
pause