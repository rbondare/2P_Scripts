@echo off
rem Run the Suite3D axon pipeline in the "suite3d" conda env.
rem   run_suite3d.bat <recording_folder> [options]      e.g.  run_suite3d.bat "Z:\Group Members\Rima\DATA_2P\AnimalRB25_261104_1200"
rem   run_suite3d.bat <recording_folder> --dry-run      only read the headers and print the plan
rem Options: see the docstring of run_suite3d_RB_axons_V1.py (or add --help).
call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate suite3d
if errorlevel 1 (
    echo Could not activate the suite3d env. Create it first: see environment_suite3d.yml
    exit /b 1
)
python "%~dp0run_suite3d_RB_axons_V1.py" %*
