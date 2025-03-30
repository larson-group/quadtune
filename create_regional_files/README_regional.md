#Quick Start

#Main Program Configuration:
The main program is CLUBB_PROF_20150117.py, where users can configure the following parameters:

Diagnostic result naming: case='custom_name'  
Output directory: outdir=path0 (for storing diagnostic results/plots)  
Simulation configuration:  
Case names: cases=[name1,name2]  
Raw data paths: filepath=[path1,path2]  
Start years: years=[2005,1979]  
Run duration: nyear=[10,1] (in years)
Model file identifier: affl=[cam,eam] (for recognizing history file naming conventions, e.g., EAM/CAM)  
History file suffixes: suffix=[h0a,h0,h1] (specifying storage parameters like h0/h1/h0a)  
Climate file paths:  
Vertical profile data: climopath=path  
Regridded data: regridpath=path  

#Plot Settings:
Resolution adjustment: pixel = 100 (modifies plot resolution)
Output format: ptype = png/pdf/eps 
Time dimension: csenson = "ANN" (for annual mean calculation)

#Subroutine Functions:

##Climate Mean Calculation
When calmean=TRUE is set, function_cal_mean.py will be called to compute climatological means.
Target Site Selection
Set latitude/longitude in main program: lats= and lons=, then enable findout=True to call function_pick_out.py  
Supports DYNCORE grids (e.g., SE/FV):  
Declare calfvsite = [True, False] (set True for FV grid)

##Search parameter: area=1.5 (matches model grids within 1.5-degree square around site, averaging multiple grids)
2D Plot Generation and Regional File Creation  
2D plots: Set draw2d=True in main program to call draw_plots_hoz_2D.py

##Regional files:
Set MKREG = True and specify resolution intll=20 or 30 in draw_hoz_plots_2D.py  
Regional files are named as: ./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc'  
Target site mean profiles will also be stored in regional files

Note: Current version can only generate complete regional files for SE grids
