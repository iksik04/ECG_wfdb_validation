@echo off
setlocal
set SCRIPT_NAME=wqrs
set DB_PATH=C:\Users\Public\Work\QRS-detector\ECG_DB\AHADB

echo Searching for records in folder %DB_PATH%...
echo.
cd %DB_PATH%
set WFDB=.

set COUNT=0

for %%f in ("%DB_PATH%\*.dat") do (
    echo Found record: %%~nf
    echo Running: %SCRIPT_NAME% -r %%~nf
    %SCRIPT_NAME% -r %%~nf
    echo.
    set /a COUNT+=1
)

echo.
echo Done! Processed records: %COUNT%
pause