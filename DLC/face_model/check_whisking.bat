@echo off
rem Blind check of the whisking on/off call (see check_whisking.py). After watching the clips and filling in
rem your_answers.csv:
rem   check_whisking.bat score --out-dir "C:\Users\rbondarenko\DLC\whisking_check\clips"

if "%~1"=="" (
    echo Usage: check_whisking.bat score --out-dir "folder with the clips"
    echo        check_whisking.bat make --out-dir "folder" "video.mp4|session_whisker.csv" [...]
    exit /b 1
)

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-tf
if errorlevel 1 exit /b 1

python -u "%~dp0check_whisking.py" %*
