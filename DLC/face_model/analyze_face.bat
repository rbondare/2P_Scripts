@echo off
rem Run the face model (nose tip + one whisker) on a raw camera video. Usage:
rem   analyze_face.bat "video.mp4" --out-dir "folder"

if "%~1"=="" (
    echo Usage: analyze_face.bat "video.mp4" --out-dir "folder"
    exit /b 1
)

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-tf
if errorlevel 1 exit /b 1

python -u "%~dp0analyze_face.py" %*
