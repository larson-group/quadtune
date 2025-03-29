# QuadTune
A tuner for global atmospheric models with a quadratic emulator.

This repo contains python source code in two folders: 
1) `create_regional_files` and 
2) `tuning_files`.

The scripts in `create_regional_files` take output from 
global atmospheric simulations and create regional files.
See `create_regional_files/README_regional.md` for more information.

The scripts in `tuning_files` take the regional files,
optimize parameter values, and display diagnostic plots.
See `tuning_files/README_tuning.md` for more information.