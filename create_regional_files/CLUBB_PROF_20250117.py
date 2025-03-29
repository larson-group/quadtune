'''
E3SM CLUBB Diagnostics package 
# python CLUBB_PROF_20250117.py
Main code to make 1) 2D plots,2) profiles, 3) budgets on selected stations, 
         and then build up  webpages  etc
    Zhun Guo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
    Kate Thayer-Calder
    Benjamin A. Stephens:stepheba@ucar.edu

'''

## ==========================================================
# Begin User Defined Settings
# User defined name used for this comparison, this will be the name 
#   given the directory for these diagnostics
case='test_202503' # A general case name
outdir='/lcrc/group/acme/public_html/diagnostic_output/ac.zguo/'

#location of the model output data
filepath=[\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\

["/data/zhun/csmruns/","/run/"],\
          ]

#location for climatology files
climopath=[\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/"],\
]

#location for regridded data
regridpath=[\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/climo/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/climo/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/climo/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/climo/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/climo/"],\
]

#location of the "run" directory with the atm_in file
runfilepath=[\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
["/lcrc/group/acme/ac.zguo/E3SM_simulations/","/run/"],\
]

cases=[ \
"sens1022_1",\
"sens1022_2",\
"sens1022_3",\
"sens1022_4",\
"sens1022_5",\
#"sens1022_6",\
#"sens1022_7",\
#"sens1022_8",\
#"sens1022_9",\
#"sens1022_10",\
#"sens1022_11",\
#"sens1022_12",\
#"sens1022_13",\
#"sens1022_14",\
#"sens1022_15",\
#"sens1022_16",\
#"sens1022_17",\
#"sens1022_18",\
#"sens1022_19",\
#"sens1022_20",\
#"sens1022_21",\
]


# Give a short name for your experiment which will appears on plots

casenames=cases#[\
#"default",\
#"Handtune",\
#"Tuner_20RG",\
#"sens1022_69",\
#"Tuner_20_",\
#"sens1022_70",\
#]
years=[\
        2005,2005, 2005,2005,2005, 2005,2005,2005, 2005,2005,2005, 2005,2005,2005, 2005,2005,2005, 2005,2005,2005, 2005,2005,2005, 2005, 1979]
nyear=[\
         5,5,5, 2,2,2, 2,2,2, 2,2,2, 2, 2,2,2, 2,2,2, 2,2,2, 2,2,2, 2]

# Affiliations, Please confirm the naming of the schema history file. Sometimes it's eam, sometimes it's cam or something else that can be declared here.
affl=['eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','eam','cam']
# suffix, Usually eam and cam use the same suffix, e.g. h0, h1, etc. But sometimes it's h0a, h1a, etc.It's better to make sure, otherwise the cal_mean function doesn't work.
suffix =['h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0','h0a']

# Although this is not required for CLUBB diagnostic package, in order to call the other diagnostic packages (E3SM or AMWG), we need to diff the SE grid to lat-lon coordinates. This is where mapfiles are needed. 
mapfile = [\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\
'/lcrc/group/acme/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc',\

'/data/zhun/inputdata/lnd/clm2/mappingdata/maps/ne30pg3/map_ne30pg3_to_0.5x0.5_nomask_aave_da_c180515.nc',\
'/data/zhun/inputdata/lnd/clm2/mappingdata/maps/ne30pg3/map_ne30pg3_to_0.5x0.5_nomask_aave_da_c180515.nc',\
]

ncopath ='/blues/gpfs/home/software/anaconda3/2021.11/envs/nco/bin/'

# NOTE, dpsc,deep scheme, has to be 'none', if silhs is turned on. 
dpsc=[\
      'zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','zm','none']

# mpsc, microphy scheme, can be P3 or MG
mpsc=[\
      'P3','P3','P3','P3','P3']

#tuning parameters:
list_of_parameters=['clubb_c7','clubb_c11','clubb_gamma_coef','clubb_c8','clubb_c_k10','clubb_c_invrs_tau_N2','clubb_c_invrs_tau_wpxp_n2_thresh','clubb_altitude_threshold','clubb_c_invrs_tau_bkgnd','clubb_c_invrs_tau_sfc','clubb_c_invrs_tau_n2_wp2','clubb_C_invrs_tau_wpxp_Ri','clubb_c_invrs_tau_n2_wpxp','clubb_c_invrs_tau_shear','clubb_c_invrs_tau_n2_xp2','clubb_c_k8','clubb_nu1','clubb_nu2','clubb_c_invrs_tau_n2_clear_wp3','clubb_c_wp2_splat','cldfrc_dp2','clubb_c_invrs_tau_n2','clubb_z_displace','micro_mg_dcs','micro_mg_autocon_lwp_exp']

# Observation Data
filepathobs='/lcrc/group/acme/ac.zguo/climatology'
#------------------------------------------------------------------------
# Setting of plots.
ptype         ='png'   # eps, pdf, ps, png, x11, ... are supported by this package
pixel         = 100    # This parameter only determines the resolution of PNG files, 
cseason       ='ANN'   # Seasons, or others
casename      =case+'_'+cseason

#------------------------------------------------------------------------
calmean          = False # make mean states
findout          = True      # pick out the locations of your sites
draw2d           = True       # 2D plots, SWCF etc.

area  = 1.5

# Note, this para helps to find out the 'ncol' within
# lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area
#------------------------------------------------------------------------
# Please give the lat and lon of sites here.
# sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   23   24   25   26   27   28   29   30   31  32   33  34   35  36
lats = [  20,  27, -20, -20,  -5,  -1,  60,   2,   9,   56,  45,   0,  10,  20,   0,   5,   9, -60,   0,   0, -45, -75,  30,  28 , 70 , 0,  -15, -15, -35, -25, -15,-15, -1,  45,  9,35]

lons = [ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311, 180, 295,  90, 205, 325, 280, 170, 340, 305,  25,  90,  90,  90, 108 , 90, 160, 260, 345, 220,  75, 180,270,270, 260,281,280]

REG  = ['DYCOMS','HAWAII','VOCAL','VOCAL_near','Namibia','Namibia_near','NP',  'SP',  'EP',  'WP',  'ITCZ', 'LBA',  'CAF', 'PA']
LatS = [    20.0,    10.0,  -25.0,     -35.0,    -20.0,          -30.0, 45.0, -60.0,   0.0, -10.0,     -10.0, -15.0,  -10.0,   -5.0]
LatN = [    35.0,    30.0,  -15.0,     -15.0,    -10.0,          -20.0, 60.0, -45.0,  15.0,  10.0,      10.0,   5.0,   10.0,   10.0]
LonW = [   226.0,   200.0,  275.0,     282.0,     -5.0,            5.0,  1.0,   1.0, 180.0, 110.0,     110.0, 285.0,   10.0,  270.0]
LonE = [   241.0,   220.0,  285.0,     288.0,      5.0,           15.0,360.0, 360.0, 260.0, 165.0,     285.0, 320.0,   40.0,  290.0]

#========================================================================

#------------------------------------------------------------------------
# Do not need to change
#------------------------------------------------------------------------

if 'P3' in mpsc and 'MG' in mpsc:
    print("NOTE:::Since both P3 and MG exist simultaneously in the microphysics, we are unable to compare them on a single plot. Therefore, drawp3 and drawmicro are forcibly set to false. ")
    drawmicrobgt = False
    drawp3 = False
if all(item == 'P3' for item in mpsc):
    drawmicrobgt = False
elif all(item == 'MG' for item in mpsc):
    drawp3 = False

ncases =len(cases)
nsite  =len(lats)

casedir=outdir+casename
print(casedir)

import numpy as np
import pdb
import os
import function_cal_mean
import function_pick_out
import draw_plots_hoz_2D

casedir=outdir+casename

if not os.path.exists(casedir):
    os.mkdir(casedir)

if calmean:
    print('Getting climatological mean')
    function_cal_mean.cal_mean(ncases, cases, years,nyear, nsite, lats, lons, area, filepath,climopath,regridpath,affl,suffix,mapfile,ncopath)

if findout:
    print('Find out the sites')
    function_pick_out.pick_out(ncases, cases, years, nsite, lats, lons, area, climopath, runfilepath, casedir,ncopath)

if draw2d:
    print('Drawing 2d')
    plot2d=draw_plots_hoz_2D.draw_2D_plot(ptype,pixel,cseason, ncases, cases, casenames, nsite, lats, lons, climopath, regridpath, runfilepath, filepathobs,casedir,list_of_parameters,REG,LatS,LatN,LonE, LonW, ncopath)
