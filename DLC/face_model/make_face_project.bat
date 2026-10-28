@echo off
rem Create the face DLC project (nose tip + one whisker) and pick the frames to label.
rem Takes ~10 min (it reads random frames from every new-setup video on the drive). Just run:
rem   make_face_project.bat

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-tf
if errorlevel 1 exit /b 1

python -u "%~dp0make_face_project.py" %*
