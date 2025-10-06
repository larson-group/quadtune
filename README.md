# QuadTune
A tuner for global atmospheric models with a quadratic emulator.

This repo contains python source code in three folders: 
1) `create_regional_files`,
2) `tuning_files` and
3) `tests`.

The scripts in `create_regional_files` take output from 
global atmospheric simulations and create regional files.
See `create_regional_files/README_regional.md` for a
user guide for creating regional files.

The scripts in `tuning_files` take the regional files,
optimize parameter values, and display diagnostic plots
on a dashboard created with plotly dash.
See `tuning_files/README_tuning.md` for a user guide for tuning.

The scripts in `tests` are pytest-tests which test the functionality of quadtune.  
To run all tests, go to the `quadtune` directory and execute `pytest`.
To run specifc tests, also go to the `quadtune` directory and execute `pytest -k ' {pattern found within one or more testfile names e.g. "v1" } '