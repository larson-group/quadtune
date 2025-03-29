This folder contains scripts that, given regional files as input,
optimize parameter values and create diagnostic plots.

Use `requirements.txt` to import needed python libraries.

In `config.py`, specify the configuration of your tuning run, 
including the names of parameters to tune, the regional metrics to use, 
and the file names containing that information.

Then run QuadTune with `python3 quadtune_driver.py` and
view the plots at http://127.0.0.1:8050/ in your web browser.

Under the hood, `quadtune_driver.py` ingests information from `config.py`
with the help of functions in `set_up_inputs.py`.  Then it
optimizes parameter values, and calls `create_nonbootstrap_figs.py`.

`create_nonbootstrap_figs.py` generates figures.

Other `bootstrap` files generate bootstrap samples
from the regional metrics.  They are experimental for now.