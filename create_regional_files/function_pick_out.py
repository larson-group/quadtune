#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:24:40 2019
    Updates on Jan 2025
    Zhun Guo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
    Kate Thayer-Calder
    Benjamin A. Stephens:stepheba@ucar.edu

"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import os
from subprocess import call


def pick_out(ncases, cases,years, nsite,lats, lons,area, climopath,runfilepath,casedir, ncopath):
# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# climopath, climatology file path
# filepathobs, filepath for observational data

# inptrs = [ncases]

 if not os.path.exists("data"):
        os.mkdir("data")

 ncea_str=ncopath+'ncea '
 ncks_str=ncopath+'ncks '
 stre="%5.2f\n"
 for im in range(0,ncases):
     infile=climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_ANN_climo.nc'
     print(infile)
     print(im)
     inptrs = Dataset(infile,'r')       # pointer to file1
     lat=inptrs.variables['lat'][:]
     nlat=len(lat)
     lon=inptrs.variables['lon'][:] 
     nlon=len(lon)
#     idx_cols=[[0 for i in range(5)] for j in range(nsite)] 
     nh=[0 for k in range(nsite)]
#     print(idx_cols)
     cols=[0,1,2,3,4]
     sits=np.linspace(0,nsite-1,nsite)

     txtfile1=runfilepath[im][0]+cases[im]+runfilepath[im][1]+'/diff*.asc'
     txtfile2=runfilepath[im][0]+cases[im]+runfilepath[im][1]+'/log*.asc'
#     os.system('mkdir '+ casedir+'/txt/')
#     os.system('cp -f '+ txtfile1+ ' '+ casedir+'/txt/')
#     os.system('cp -f '+ txtfile2+ ' '+ casedir+'/txt/')


     os.system('rm -f ./data/'+cases[im]+'_site_location.nc')
     outf =Dataset('./data/'+cases[im]+'_site_location.nc','w')
     outf.createDimension("sit",nsite)
     outf.createDimension("col",5)
#     outf.variables['sit'][:]=sits
#     outf.variables['col'][:]=cols
     outf.createVariable('idx_cols','i',('sit','col'))
     outf.createVariable('n','i',('sit'))
     outf.variables['n'][:]=0
     outf.variables['idx_cols'][:,:]=0

# ========================================================================== 
# find out the cols and their numbers
#    the axis of site is stored in idx_cols(site,n)
#     localtxt=''
     os.system('set sw_ea = 0')
     os.system('set num = 0')
     for i in range(0,nlat):
#         if (lon[i] >= 230) & (lon[i] < 241) & (lat[i]>=23) & (lat[i] < 33):
#             localtxt=localtxt+','+str(i)
#             os.system(ncks_str+'-d ncol,'+str(i)+' '+infile +' -O '+'./data/tmpcol_'+str(i)+'.nc')
#              os.system('set sw = `ncks -H -s "%5.2f\n"  -d ncol,'+str(i)+' -v SWCF '+infile + '`')
#              os.system('echo $sw') 
#              os.system(  'set \!:1 = `echo  | bc -l` sw_ea = $sw_ea + $sw')
##              os.system(  'set \!:1 = `echo  | bc -l` num = 1 + $num')
         for j in range(0,nsite): 
             if (lon[i] >= lons[j]-area) & (lon[i] < lons[j]+area) & (lat[i]>=lats[j]-area) & (lat[i] < lats[j]+area): 
                 outf.variables['idx_cols'][j,nh[j]]=i 
                 outf.variables['n'][j]=nh[j]+1
#     os.system(ncks_str+'-d ncol'+localtxt+' '+infile +' -O '+'./data/sc.nc')
#     os.system(ncea_str+' ./data/tmpcol_*.nc' + '-O '+'./data/sc.nc')
#     os.system('rm -f ./data/tmpcol_*.nc' + '-O '+'./data/sc.nc')
     outf.close()
