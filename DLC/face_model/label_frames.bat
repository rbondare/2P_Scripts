@echo off
rem Open the face frames (nose tip + one whisker) in the labelling GUI. Just run:   label_frames.bat

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-gui
if errorlevel 1 exit /b 1

python -u "%~dp0label_frames.py" %*
