# User guide for creation of regional files

This folder contains scripts that, given global-model output,
create the regional files that are later fed into the tuning scripts.
The regional files are assumed to be in netcdf format,
with one regional file per each of the 2P+1 global model runs.
A regional file contains all of the regional metrics ---
e.g., shortwave cloud forcing from stratocumulus region 6_14 ---
from a single run.  

## Quick Start

### Main Program Configuration:

The main program is CLUBB_PROF_20150117.py.  In this file, users 
can configure the following parameters:

Diagnostic result naming: case='custom_name'  
Output directory: outdir=path0 (for storing diagnostic results/plots)  
Simulation configuration:  
- Names of default and sensitivity simulations: cases=[name1,name2]  
- Path of the raw model output: filepath=[path1,path2]  
- Start years of simulations: years=[2005,1979]  
- Run durations of simulations: nyear=[10,2] (in years, 14-month minimum)
- Model file identifier: affl=[cam,eam] (for recognizing history file naming conventions, e.g., EAM/CAM)  
- History file suffixes: suffix=[h0a,h0,h1] (specifying history storage stream)  

Model output climatology file paths:  
- Vertical profile data: climopath=path  
- Regridded output data: regridpath=path  

### Plot Settings:

Resolution adjustment: pixel = 100 (modifies plot resolution)  
Output format: ptype = png/pdf/eps  
Time dimension: cseason = "ANN" (for annual mean calculation)

## Subroutine Functions:

### Climate Mean Calculation

When calmean=TRUE is set, function_cal_mean.py will be called to compute climatological means.  
Target Site (Custom regions) Selection:  
- Set latitude/longitude in main program: lats= and lons=, then enable findout=True to call function_pick_out.py  


### Search parameter for location of custom regions: 

Proximity of region to grid column: area = 1.5 (matches model grids within 1.5-degree square around site, averaging multiple grids)  
2D Plot Generation and Regional File Creation:  
- 2D plots: Set draw2d=True in main program to call draw_plots_hoz_2D.py  

### Regional files:

To create regional files:
- Set MKREG = True and 
- Specify resolution (regional box size) intll=20 or 30 in draw_hoz_plots_2D.py  
Regional files are named as: ./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc'  
Target site mean profiles will also be stored in regional files  

Note: Current version can only generate complete regional files for SE grids, not FV grids.
