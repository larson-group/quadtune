#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Updated on Jan 17 2025

@authors: 
    Zhun Guo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
    Kate Thayer-Calder
    Benjamin A. Stephens:stepheba@ucar.edu
"""

from netCDF4 import Dataset
import numpy as np
import scipy as sp
import os

from subprocess import call

def cal_mean(ncases, cases,years,nyear, nsite,lats, lons,area, filepath, climopath, regridpath,affl, suffix, mapfile, ncopath):
# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data

 ncea_str    =ncopath+'ncea '
 ncwa_str    =ncopath+'ncwa -C '
 ncks_str    =ncopath+'ncks '
 ncap2_str   =ncopath+'ncap2 '
 ncrcat_str  =ncopath+'ncrcat '
 ncrename_str=ncopath+'ncrename '

 for im in range (0, ncases ):
#  call('rm -f  ',climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_*_climo.nc')
     datalocal = filepath[im][0] +cases[im]+filepath[im][1]
     print(datalocal)
     os.system('mkdir -p '+climopath[im][0]+cases[im]+climopath[im][1])
     os.system('mkdir -p '+regridpath[im][0]+cases[im]+regridpath[im][1])
   
     outfile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_DJF_climo.nc'
     infile=' '
     for yr in range (0, nyear[im] ):
         infile=infile+datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1+yr).rjust(4,'0')+'-01.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1+yr).rjust(4,'0')+'-02.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-12.nc '

     os.system(ncea_str+infile +' -O '+outfile)
     climofile=regridpath[im][0]+cases[im]+regridpath[im][1]+cases[im]+'_DJF_climo.nc'
     
     
     os.system(ncks_str +' --map='+mapfile[im]+'   ' + outfile +' -O ' +climofile)   

     outfile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_MAM_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-03.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-04.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-05.nc '
     os.system(ncea_str+infile +' -O '+outfile)
   
     outfile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_JJA_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-06.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-07.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-08.nc '

     os.system(ncea_str+infile +' -O '+outfile)

     climofile=regridpath[im][0]+cases[im]+regridpath[im][1]+cases[im]+'_JJA_climo.nc'
     
     os.system(ncks_str +' --map='+mapfile[im]+'   ' + outfile +' -O ' +climofile)  

     outfile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_SON_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-09.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-10.nc '+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+yr).rjust(4,'0')+'-11.nc '

     os.system('rm -f '+climopath[im][0]+cases[im]+climopath[im][1]+'*_ANN_*.nc')
     os.system(ncea_str+infile +' -O '+outfile)
   
     outfile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_climo.nc'
     climofile=regridpath[im][0]+cases[im]+regridpath[im][1]+cases[im]+'_ANN_climo.nc'
     infile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_SON_climo.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_JJA_climo.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_MAM_climo.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_DJF_climo.nc '
     os.system(ncea_str+infile +' -O '+outfile)
     

     os.system(ncks_str +' --map='+mapfile[im]+'   ' + outfile +' -O ' +climofile) 


#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,-5,5 -d lon,290,300 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_LBA_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_LBA_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_LBA.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_LBA -v .SWCF,SWCF_LBA -v .LWCF,LWCF_LBA -v .CLDTOT,CLDTOT_LBA -v .TMQ,TMQ_LBA '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_LBA.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,23,33 -d lon,230,241 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_DYCOMS_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_DYCOMS_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_DYCOMS.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_DYCOMS -v .SWCF,SWCF_DYCOMS -v .LWCF,LWCF_DYCOMS -v .CLDTOT,CLDTOT_DYCOMS -v .TMQ,TMQ_DYCOMS '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_DYCOMS.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,15,25 -d lon,200,210 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_HAWAII_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_HAWAII_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_HAWAII.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_HAWAII -v .SWCF,SWCF_HAWAII -v .LWCF,LWCF_HAWAII -v .CLDTOT,CLDTOT_HAWAII -v .TMQ,TMQ_HAWAII '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_HAWAII.nc')
#
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,-25,-15 -d lon,275,280 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_VOCAL_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_VOCAL_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_VOCAL.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_VOCAL -v .SWCF,SWCF_VOCAL -v .LWCF,LWCF_VOCAL -v .CLDTOT,CLDTOT_VOCAL -v .TMQ,TMQ_VOCAL '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_VOCAL.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,-5,5 -d lon,110,140 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_WP_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_WP_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_WP.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_WP -v .SWCF,SWCF_WP -v .LWCF,LWCF_WP -v .CLDTOT,CLDTOT_WP -v .TMQ,TMQ_WP '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_WP.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,0,15 -d lon,180,270 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_EP_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_EP_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_EP.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_EP -v .SWCF,SWCF_EP -v .LWCF,LWCF_EP -v .CLDTOT,CLDTOT_EP -v .TMQ,TMQ_EP '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_EP.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,-60,-45 -d lon,180,200 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SP_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SP_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SP.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_SP -v .SWCF,SWCF_SP -v .LWCF,LWCF_SP -v .CLDTOT,CLDTOT_SP -v .TMQ,TMQ_SP '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SP.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,45,60 -d lon,180,200 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_NP_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_NP_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_NP.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_NP -v .SWCF,SWCF_NP -v .LWCF,LWCF_NP -v .CLDTOT,CLDTOT_NP -v .TMQ,TMQ_NP '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_NP.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,0,5 -d lon,270,280 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_PA_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_PA_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_PA.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_PA -v .SWCF,SWCF_PA -v .LWCF,LWCF_PA -v .CLDTOT,CLDTOT_PA -v .TMQ,TMQ_PA '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_PA.nc')
#
#     os.system(ncwa_str+' -O -C -a lat,lon,time  -d lat,-5,5 -d lon,20,30 '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SAF_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SAF_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SAF.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_SAF -v .SWCF,SWCF_SAF -v .LWCF,LWCF_SAF -v .CLDTOT,CLDTOT_SAF -v .TMQ,TMQ_SAF '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SAF.nc')
#
#     os.system(ncwa_str+'  --flt -O -a lat,lon -w area  '+ climofile + ' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB_ALL.nc')
#     os.system(ncks_str+' -v PRECL,SWCF,LWCF,CLDTOT,TMQ '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB_ALL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc')
#     os.system(ncrename_str+' -v .PRECL,PRECT_GLB -v .SWCF,SWCF_GLB -v .LWCF,LWCF_GLB -v .CLDTOT,CLDTOT_GLB -v .TMQ,TMQ_GLB '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc')
#
#
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SAF.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_PA.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_WP.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_NP.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_SP.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_HAWAII.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_DYCOMS.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_VOCAL.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_LBA.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#
#     os.system(ncks_str+'-A '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_EP.nc '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc' )
#
#     with open(climopath[im][0]+cases[im]+'/run/atm_in') as inf:
#         clubb_c_invrs_tau_n2 = 0
#         for line in inf:
#             line = line.split('=')
#             line[0] = line[0].strip()
#             if line[0] == 'clubb_c_invrs_tau_n2':
#                 clubb_c_invrs_tau_n2 = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_n2_clear_wp3':
#                 clubb_c_invrs_tau_n2_clear_wp3 = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_n2_wpxp':
#                 clubb_c_invrs_tau_n2_wpxp = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_n2_xp2':
#                 clubb_c_invrs_tau_n2_xp2 = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_n2_wp2':
#                 clubb_c_invrs_tau_n2_wp2 = float(line[1])
#             elif line[0] == 'clubb_c_k10h':
#                 clubb_c_k10h = float(line[1])
#             elif line[0] == 'clubb_c_k10':
#                 clubb_c_k10 = float(line[1])
#             elif line[0] == 'clubb_c8':
#                 clubb_c8 = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_wpxp_ri':
#                 clubb_c_invrs_tau_wpxp_ri = float(line[1])
#             elif line[0] == 'clubb_c_invrs_tau_wpxp_n2_thresh':
#                 clubb_c_invrs_tau_wpxp_n2_thresh = float(line[1])
#     os.system(ncap2_str +'-s clubb_c_invrs_tau_n2='+str(clubb_c_invrs_tau_n2)+';clubb_c_invrs_tau_n2_clear_wp3='+str(clubb_c_invrs_tau_n2_clear_wp3)+';clubb_c_invrs_tau_n2_wpxp='+str(clubb_c_invrs_tau_n2_wpxp)+';clubb_c_invrs_tau_n2_xp2='+str(clubb_c_invrs_tau_n2_xp2)+';clubb_c_invrs_tau_n2_wp2='+str(clubb_c_invrs_tau_n2_wp2)+';clubb_c_k10h='+str(clubb_c_k10h)+';clubb_c_k10='+str(clubb_c_k10)+';clubb_c8='+ str(clubb_c8)+';clubb_c_invrs_tau_wpxp_ri='+str(clubb_c_invrs_tau_wpxp_ri)+'clubb_c_invrs_tau_wpxp_n2_thresh='+str(clubb_c_invrs_tau_wpxp_n2_thresh) +' '+climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_GLB.nc '+' ./data/'+cases[im]+'_ANN_Regionalmean.nc')
     ln_str='ln -s '
     os.system(ln_str+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_SON_climo.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_SON_budget_climo.nc' )
     os.system(ln_str+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_MAM_climo.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_MAM_budget_climo.nc' )
     os.system(ln_str+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_JJA_climo.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_JJA_budget_climo.nc' )
     os.system(ln_str+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_DJF_climo.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_DJF_budget_climo.nc' )
     os.system(ln_str+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_climo.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-10.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_10_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-10.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_10_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-11.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_11_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-11.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_11_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-12.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_12_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-12.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_12_budget_climo.nc' )
   
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1).rjust(4,'0')+'-01.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_01_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1).rjust(4,'0')+'-01.nc'  +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_01_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1).rjust(4,'0')+'-02.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_02_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]+1).rjust(4,'0')+'-02.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_02_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-03.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_03_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-03.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_03_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-04.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_04_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-04.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_04_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-05.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_05_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-05.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_05_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-06.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_06_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-06.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_06_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-07.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_07_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-07.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_07_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-08.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_08_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-08.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_08_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-09.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_09_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.'+affl[im]+'.'+suffix[im]+'.'+str(years[im]).rjust(4,'0')+'-09.nc' +' '+ climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_09_budget_climo.nc' )
