@echo off
rem Run the body + eye DeepLabCut pipeline on recordings stored on the drive:
rem copies the raw video to a local work folder, runs DLC there, writes the quality-metric report, then
rem copies body_results\, eye_results\ and success.log back next to the raw video.
rem It also writes whisking on/off, amplitude and face motion to whisker_results\ (face_model\face_motion.py);
rem add --no-whisker to skip that.
rem See run_dlc_to_drive.py for details.
rem
rem   run_dlc.bat "Z:\Group Members\Rima\Camera\AnimalRB22\AnimalRB22_261023_1427"
rem   run_dlc.bat --scan "Z:\Group Members\Rima\Camera" --dry-run
rem   run_dlc.bat --scan "Z:\Group Members\Rima\Camera"
rem
rem The DLC models are not in the repo: experiment_config.yaml points to C:\Users\rbondarenko\DLC\Training_Data.

if "%~1"=="" (
    echo Usage: run_dlc.bat "session_folder" [more session folders]
    echo        run_dlc.bat --scan "root_folder" [--dry-run]
    exit /b 1
)

call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate dlc-tf
if errorlevel 1 exit /b 1

python -u "%~dp0run_dlc_to_drive.py" %*
