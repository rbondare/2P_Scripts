@echo off
rem Train + evaluate the face model (nose tip + one whisker) on the hand labels. Just run:   train_face_model.bat
rem Options are passed on to train_face_model.py, e.g.   train_face_model.bat --maxiters 60000

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-tf
if errorlevel 1 exit /b 1

python -u "%~dp0train_face_model.py" %*
