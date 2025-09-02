# User guide for tuning

This folder contains scripts that, given regional files as input,
optimize parameter values and create diagnostic plots. 
The scripts have been tested using Python 3.12.
The regional files must be created beforehand by running the scripts
in folder `create_regional_files`.

## Getting started

1) Use `requirements.txt` to import needed python libraries.
E.g., use your IDE or at the command line type 
`pip install -r requirements.txt`.

2) QuadTune has no namelist file, but in `config_default.py`, 
you may specify the configuration of your tuning run, 
including the names of parameters to tune, the regional metrics to use, 
and the file names containing that information.
For more information on the configurable quantities,
see the code comments in `config_default.py`.
The file `config_example.py` is a scratch file
where a user can store old configurations.

3) Then run QuadTune with `python3 quadtune_driver.py` and
view the plots at http://127.0.0.1:8050/ in your web browser.

## Some setup variables in config_default.py

`folder_name`:  This designates the folder where the regional files are stored.  There are 2P regional .nc files plus a default file.
`obsPathAndFilename`: This designates the folder and name of the .nc file where observational data is recorded.

## Under the hood

QuadTune takes  information from regional files in netcdf format.

Under the hood, `quadtune_driver.py` ingests information from `config.py`
with the help of functions in `set_up_inputs.py`.  Then it
optimizes parameter values, and calls `create_nonbootstrap_figs.py`.

`create_nonbootstrap_figs.py` generates figures and displays
figures on a plotly dash dashboard.

The other files, namely `create_bootstrap_figs.py` and `do_bootstrap_calcs.py`,
generate bootstrap samples from the regional metrics and create plots.
These scripts are experimental for now.
