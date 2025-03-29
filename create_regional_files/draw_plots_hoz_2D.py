'''
    Draw 2D plots.Set MKREG=True, to create Regional files (set Intll=xx), the profiles of selected sites are also contained in RG files.
    Updates on Jan 2025 
    Zhun Guo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
    Kate Thayer-Calder
    Benjamin A. Stephens:stepheba@ucar.edu
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import math
import os
import Common_functions
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cmaps

from matplotlib import font_manager as fm
from subprocess import call
from scipy.interpolate import griddata

def convert_longitude(longitude):
#    if longitude > 180:
        return longitude - 180.
#    else:
#        return longitude 

def traditional_round(number, ndigits=0):
    multiplier = 10 ** ndigits
    return int((number * multiplier + 0.5) / multiplier)

def mse(obs,pre):
    obs=obs.flatten()
    pre=pre.flatten()
    return np.sum((obs-pre)**2)/len(obs)

def rmsq(obs,pre):
    obs=obs.flatten()
    pre=pre.flatten()
    return (obs-pre)**2/len(obs)

def get_varible_name(var_org):
    for item in sys._getframe().f_locals.items():
        print(item[0],item[1])
    for item in sys._getframe(1).f_locals.items():
        print(item[0],item[1])
    for item in sys._getframe(2).f_locals.items():
        if (var_org is item[1]):
            return item[0]

def get_name(number):
    print("{} = {}".format(get_varible_name(number),number))

def rmse(obs,pre):
    obs=obs.flatten()
    pre=pre.flatten()
    return np.sqrt(np.sum((obs-pre)**2)/len(obs))

def rmsep(obs,pre):
    obs=obs.flatten()
    pre=pre.flatten()
    obsmean=np.mean(obs)
    premean=np.mean(pre)
    obs=obs-obsmean
    pre=pre-premean
    return np.sqrt(np.sum((obs-pre)**2)/len(obs))

def racc(obs,pre):
    obs=obs.flatten()
    pre=pre.flatten()
    obsmean=np.mean(obs)
    premean=np.mean(pre)
    obs=obs-obsmean
    pre=pre-premean
    return 1 - np.sum(obs*pre)/len(obs)/np.sqrt(np.sum((obs)**2)/len(obs))/np.sqrt(np.sum((pre)**2)/len(obs))  #/len(obs)

def draw_2D_plot (ptype,pixel,cseason, ncases, cases, casenames, nsite, lats, lons, climopath, regridpath, runfilepath, filepathobs,casedir,list_of_parameters,REG,LatS,LatN,LonE, LonW, ncopath):

# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data
# inptrs = [ncases]
 if not os.path.exists(casedir):
        os.mkdir(casedir)

 if not os.path.exists(casedir+'/2D'):
        os.mkdir(casedir+'/2D') 

 intll   = 30
 MKREG   = True
 edgelon = 360/int(intll) + 1
 edgelat = 180/int(intll) + 1
 clevel  = 500

 infiles  = ['' for x in range(ncases)] 
 infiles2 = ['' for x in range(ncases)]
 infiles3 = ['' for x in range(ncases)]
 ncdfs    = ['' for x in range(ncases)] 
 alpha    = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','A','B','C','D','E','F','G','H','I','J','K','L','M','N']
 cunits = ['']

 varis    = ['SWCF',          'MSWCF',         'LWCF',           'PRECT','O500' ,'LHFLX','SHFLX',  'TMQ','PSL','TS', 'U10', 'CLDTOT'    , 'CLDLOW'   ,    'CLDHGH']
 varisobs = ['toa_cre_sw_mon','toa_cre_sw_mon','toa_cre_lw_mon', 'PRECT','wap'  ,'hfls', 'hfss',   'prw','psl','tas','uas', 'CLDTOT_CAL','CLDLOW_CAL','CLDHGH_CAL'] #,'hfls','hfss','hfls','hfss']
 cscale   = [     1,                1,               1,         86400000,      1,     1,       1,      1,    1,    1,    1,     100,         100,         100,       1000,1,1,1,1,1,1,1]
 cscaleobs =  [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
 varisl =  ['SWCF','MSWCF','LWCF', 'PRECT','O500','RESTOM','TMQ','PSL','CLDTOT','U10','LHFLX','SHFLX'] #,'wprtp_sfc','wpthlp_sfc','upwp_sfc','vpwp_sfc']
 obsdataset=['ceres_ebaf_toa_v4.1','ceres_ebaf_toa_v4.1','ceres_ebaf_toa_v4.1','GPCP_v2.3',  'ERA-Interim',   'ERA-Interim',  'ERA-Interim', 'ERA-Interim',  'ERA-Interim',  'ERA-Interim','ERA-Interim','CALIPSOCOSP', 'CALIPSOCOSP','CALIPSOCOSP','ERA-Interim',    'ERA-Interim',  'ERA-Interim',  'ERA-Interim']
 rangeyr   =[      '200101_201812',      '200101_201812',     '200101_201812',          '','197901_201612', '197901_201612','197901_201612','197901_201612','197901_201612','197901_201612','197901_201612','','','','197901_201612','197901_201612','197901_201612','197901_201612']

# no need to edit
 ncea_str    = ncopath+'/ncea '
 ncwa_str    = ncopath+'/ncwa -C '
 ncks_str    = ncopath+'/ncks '
 ncap2_str   = ncopath+'/ncap2 '
 ncecat_str  = ncopath+'/ncecat '
 ncrename_str= ncopath+'/ncrename '
 ncdiff_str  = ncopath+'/ncdiff '

 nvaris = len(varis)
 plot2d=['' for x in range(nvaris)]
 cmap1 = cmaps.MPL_rainbow 

 for iv in range(0, nvaris):
# make plot for each field 
   if(varis[iv] == 'CLDTOT' or varis[iv] == 'CLDLOW' or varis[iv] == 'CLDHGH'):
       levels = [ 2, 5, 10, 20, 30, 40, 50, 60, 70,80, 90]
   if(varis[iv] == 'LWCF'):
       levels = [0, 10,  20, 30,  40, 50, 60, 70, 80, 90] 
   if(varis[iv] == 'SWCF' or varis[iv] =='FLUT'):
       levels = [ -130,-120, -110, -100, -90, -80, -70, -60, -50,-40,-30,-20,-10,0]
   if(varis[iv] == 'PRECT' or varis[iv]=='QFLX'):
       levels = [0.5,1.0,2,  3, 4,5, 6,7,8, 9, 10,  12,13,14,15,16]
   if(varis[iv] == 'LHFLX'):
       levels = [0,5, 15, 30, 60, 90, 120, 150, 180,210, 240, 270, 300]
   if(varis[iv] == 'SHFLX'):
       levels = [-100,-75, -50, -25, -10, 0, 10, 25, 50, 75, 100, 125,150]
   if(varis[iv] == 'U10'):
       levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
   if(varis[iv] == 'TMQ'):
       levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
   if(varis[iv] == 'TGCLDLWP'):
       levels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]

#  Observational data
   if(obsdataset[iv] =='CCCM'):
       if(cseason == 'ANN'):
           fileobs = '/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-'+cseason+'.nc'
       else:
           fileobs = '/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-2010-'+cseason+'.nc'
   else:
       if (varisobs[iv] == 'PRECT' or varis[iv] == 'CLDTOT' or varis[iv] == 'CLDLOW' or varis[iv] == 'CLDHGH'):
           fileobs = filepathobs+'/'+obsdataset[iv]+'/'+obsdataset[iv]+'_'+cseason+'_climo.nc'
       else:
           fileobs = filepathobs + '/'+obsdataset[iv]+'/'+obsdataset[iv]+'_'+cseason+'_'+rangeyr[iv] +'_climo.nc'

   inptrobs = Dataset(fileobs,'r') 
   latobs=inptrobs.variables['lat'][:]
   lonobs=inptrobs.variables['lon'][:]
   if (varisobs[iv] =='U10'):
      B0=inptrobs.variables[varisobs[iv]][0,:,:] 
      B1=inptrobs.variables['V10'][0,:,:]
      B=(B0*B0+B1*B1)
      B=B * cscaleobs[iv]
      B=np.sqrt(B)
   elif (varisobs[iv] == 'MSWCF'):
      B=inptrobs.variables['toa_cre_sw_mon'][0,:,:]
      B=B * cscaleobs[iv] *(-1)
   elif (varisobs[iv] == 'wap'):
      levobs=inptrobs.variables['plev'][:]
      levobs_idx = np.abs(levobs - clevel*100).argmin()
      B=inptrobs.variables['wap'][0,levobs_idx,:,:]
   else:
      B=inptrobs.variables[varisobs[iv]][0,:,:]
      B=B * cscaleobs[iv]
      if (obsdataset[iv] =='CALIPSOCOSP'):
          B=np.maximum(B,0)
   inptrobs.close()
# Obs data

   if MKREG:
       # Preparation of provisional datafile based on the most original observations for the calculation of global and regional averages. These files (local_*) will be deleted after completion of the calculations. 
        os.system('cp -rf ' + fileobs + ' local_obs0.nc')
        os.system(ncap2_str+' -h -O -s  "weights=cos(lat*3.1415/180)"  local_obs0.nc  local_obs.nc')
#        print(ncap2_str+' -h -O -s  "weights=cos(lat*3.1415/180)"  local_obs0.nc  local_obs.nc')
        if (varisobs[iv] == 'MSWCF'):
            os.system(ncap2_str+' -s "toa_cre_sw_mon=toa_cre_sw_mon*(-1);" local_obs0.nc  local_obs.nc')
        if (varisobs[iv] == 'O500'):
            os.system(ncks_str+' -v  '+varisobs[iv] +' local_obs0.nc  local_obs1.nc')
            os.system(ncea_str+' -F -d plev,16,16 local_obs1.nc  local_obs.nc' )
        # the Temporary file local_obs.nc should have same size with the most original observations!

        # 1 making global means
        tmp_str = ncwa_str+' -O -w weights -a lat,lon local_obs.nc  glb0_'+str(intll)+'_OBS.nc'
        os.system(tmp_str)
        
        # If this is the first variable, then a general glb file will be generated and subsequent variables will all be consolidated in this file, namly str(intll)+'_OBS.nc'.
        if  (iv == 0):
            os.system('rm -f '+ str(intll)+'_OBS.nc' )
            os.system(ncrename_str+' -v '+varisobs[iv] +','+varis[iv]+'_GLB  glb0_' +str(intll)+'_OBS.nc ' + str(intll)+'_OBS.nc' )
        else:
            os.system(ncrename_str+' -v '+varisobs[iv] +','+varis[iv]+'_GLB  glb0_' +str(intll)+'_OBS.nc glb_' + str(intll)+'_OBS.nc' )
            os.system(ncks_str +' -A glb_'+str(intll)+'_OBS.nc '+str(intll) +'_OBS.nc')
        os.system('rm -f glb*.nc')

        # 2 making specifc RG means
        for ir in range(0, len(REG)): 
            tmp_str = ncwa_str+' -O -w weights -a lat,lon  -d lat,'+str(LatS[ir])+','+str(LatN[ir])+' -d lon,'+str(LonW[ir])+','+str(LonE[ir])+' local_obs.nc srg.nc' 
            os.system(tmp_str)
            os.system(ncks_str+' -v '+varisobs[iv] + ' srg.nc -O '+' srg_'+REG[ir]+'.nc')
            os.system(ncrename_str+' -v '+varisobs[iv] +','+varis[iv]+'_'+REG[ir]+' srg_'+REG[ir]+'.nc  ' + 'local_'+REG[ir]+'.nc' )

            os.system('rm -f srg*.nc')
            os.system(ncks_str +' -A local_'+REG[ir]+'.nc '+str(intll) +'_OBS.nc')

        for slat in range(1, int(edgelat)):
            for slon in range(1, int(edgelon)):
                tmp_str = ncwa_str+' -O  -w weights -a lat,lon  -d lat,'+str(90.-slat*intll)+','+str(90.+intll-slat*intll)+' -d lon,'+str(1.0*slon*intll-intll)+','+str(1.0*slon*intll)+' local_obs.nc latlon.nc'
                os.system(tmp_str)
                os.system(ncks_str+' -v '+varisobs[iv] + ' latlon.nc -O '+' latlon_'+str(slat)+'_'+str(slon)+'.nc')
                os.system(ncrename_str+' -v '+varisobs[iv] +','+varis[iv]+'_'+str(slat)+'_'+str(slon)+' latlon_'+str(slat)+'_'+str(slon)+'.nc  ' + 'local_'+str(slat)+'_'+str(slon)+'.nc' )
                os.system(ncks_str +' -A local_'+str(slat)+'_'+str(slon)+'.nc ' + str(intll)+'_OBS.nc')

                tmp_str = ncwa_str+' -O  -a lat,lon -d lat,'+str(90.-slat*intll)+','+str(90.+intll-slat*intll)+' -d lon,'+str(1.0*slon*intll-intll)+','+str(1.0*slon*intll)+' local_obs.nc latlon2.nc'
                os.system(tmp_str)
                os.system(ncks_str+' -v weights latlon2.nc -O '+' latlon2_'+str(slat)+'_'+str(slon)+'.nc')
                os.system(ncrename_str+' -v weights,weights_'+str(slat)+'_'+str(slon)+'_'+varis[iv]+' latlon2_'+str(slat)+'_'+str(slon)+'.nc  ' + 'local_2_'+str(slat)+'_'+str(slon)+'.nc' )
                os.system(ncks_str +' -A local_2_'+str(slat)+'_'+str(slon)+'.nc ' + str(intll)+'_OBS.nc')

                os.system('rm -f latlon*.nc') 
                
        os.system('rm -f local_*.nc')
 
   #************************************************
   # create plot
   #************************************************
   plotname = casedir+'/2D/Horizontal_'+varis[iv]+'_'+cseason
   plot2d[iv] = 'Horizontal_'+varis[iv]+'_'+cseason

   nplt     = ncases+1

   fig, axes = plt.subplots( nrows= traditional_round(nplt/2), ncols=2, figsize=(12, 8 ), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
   axes      = axes.flatten()

   for im in range(0, ncases):
       ncdfs[im]  = './data/'+cases[im]+'_site_location.nc' 
       infiles[im] = climopath[im][0]+cases[im]+climopath[im][1]+cases[im]+'_'+cseason+'_climo.nc'
       inptrs  = Dataset(infiles[im],'r')       

       lat=inptrs.variables['lat'][:]
       nlat=len(lat)
       lon=inptrs.variables['lon'][:]
       nlon=len(lon)
       area=inptrs.variables['area'][:]

       area_wgt = np.zeros(nlat)

       sits=np.linspace(0,nsite-1,nsite)
       ncdf= Dataset(ncdfs[im],'r')
       n   =ncdf.variables['n'][:]
       idx_cols=ncdf.variables['idx_cols'][:]

       infiles2[im] = regridpath[im][0]+cases[im]+regridpath[im][1]+cases[im]+'_ANN_climo.nc'
       inptrs2 = Dataset(infiles2[im],'r')       # pointer to file2
       rlat=inptrs2.variables['lat'][:]
       rlon=inptrs2.variables['lon'][:]

       if (varis[iv] != 'O500' and varis[iv] != 'PRECT' and varis[iv] != 'FLUT' and varis[iv] != 'U10' and varis[iv] != 'MSWCF'):
           A  = inptrs.variables[varis[iv]][0,:]
           C  = inptrs2.variables[varis[iv]][0,:,:]
       elif (varis[iv] == 'MSWCF'):
           A  = inptrs.variables['SWCF'][0,:]
           A  = A* (-1)
           C  = inptrs2.variables['SWCF'][0,:,:]
           C  = C* (-1)
       elif (varis[iv] == 'PRECT'):
           A  = inptrs.variables['PRECC'][0,:]+inptrs.variables['PRECL'][0,:]
           C  = inptrs2.variables['PRECC'][0,:,:]+inptrs2.variables['PRECL'][0,:,:]
       elif (varis[iv] == 'U10'):
           A  = inptrs.variables['U10'][0,:]*inptrs.variables['U10'][0,:]
           A  = np.sqrt(A)
           C  = inptrs2.variables['U10'][0,:,:]*inptrs2.variables['U10'][0,:,:]
           C  = np.sqrt(C)
       elif (varis[iv] == 'O500'):
           lev=inptrs.variables['lev'][:]
           lev_idx = np.abs(lev - clevel).argmin()
           A  = inptrs.variables['OMEGA'][0,lev_idx,:]
       elif (varis[iv] == 'FLUT'):
           A  = inptrs.variables['FLUT'][0,:]-inptrs.variables['FLNS'][0,:]
           C  = inptrs2.variables['FLUT'][0,:,:]-inptrs2.variables['FLNS'][0,:,:]

       A_xy = A
       A_xy = A_xy * cscale[iv]
       C_xy = C
       C_xy = C_xy * cscale[iv]

       parameter_values=np.zeros(len(list_of_parameters))

       if MKREG:
          if (iv == 0):
             if (os.path.exists('./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc')):
                 os.system('/bin/rm -f '+'./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc')

             for index,parameter in enumerate(list_of_parameters):
                 with open(runfilepath[im][0]+cases[im]+runfilepath[im][1]+'atm_in') as inf:
                     for line in inf:
                         line = line.split('=')
                         line[0] = line[0].strip()
                         if line[0] == parameter:
                             #print(parameter," = ",parameter_values[index])
                             parameter_values[index]=float(line[1])

             outf = Dataset('./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc','w')
             outf.createDimension("col",1)
             for ivl in range(0,len(varisl)):
                 #print(varisl[ivl])
                 lev=inptrs.variables['lev'][:]
                 lev_idx1 = np.abs(lev - 200).argmin()
                 lev_idx2 = np.abs(lev - 500).argmin()
                 lev_idx3 = np.abs(lev - 700).argmin()
                 locals()[varisl[ivl]+'_GLB'] = 0
                 outf.createVariable(varisl[ivl]+'_GLB','f4',('col'))

                 if (varisl[ivl] == "PRECT"):
                      locals()[varisl[ivl]] = (inptrs.variables['PRECC'][0,:]+inptrs.variables['PRECL'][0,:])*24*3600*1000
                 elif(varisl[ivl] == "RESTOM"):
                      locals()[varisl[ivl]] = inptrs.variables['FLNT'][0,:]+inptrs.variables['FSNT'][0,:]
                 elif(varisl[ivl] == "MSWCF"):
                      locals()[varisl[ivl]] = inptrs.variables['SWCF'][0,:]*(-1)
                 elif(varisl[ivl] == "O200" or varisl[ivl] == "O500" or varisl[ivl] == "O700"):
                      OMEGA = inptrs.variables['OMEGA'][0,:,:]
                      O200=OMEGA[lev_idx1,:]
                      O500=OMEGA[lev_idx2,:]
                      O700=OMEGA[lev_idx3,:]
                 elif (varisl[ivl] == "U200" or varisl[ivl] == "U500" or varisl[ivl] == "U700"):
                      U = inptrs.variables['U'][0,:,:]
                      U200=U[lev_idx1,:]
                      U500=U[lev_idx2,:]
                      U700=U[lev_idx3,:]
                 else:
                      locals()[varisl[ivl]] = inptrs.variables[varisl[ivl]][0,:]

                 if (varisl[ivl] == 'SWCF'):
                      os.system(ncks_str+' -v toa_cre_sw_mon /home/ac.zguo/diagnostic_v2_0/climo/ceres_ebaf_toa_v4.1/ceres_ebaf_toa_v4.1_ANN_200101_201812_climo.nc -O '+cases[im]+'_local_obs0.nc')
                      os.system(ncrename_str+' -v toa_cre_sw_mon,SWCF '+cases[im]+'_local_obs0.nc '+cases[im]+'_local_obs1.nc')
                      os.system(ncdiff_str+regridpath[im][0]+cases[im]+regridpath[im][1]+cases[im]+'_ANN_climo.nc '+cases[im]+'_local_obs1.nc '+cases[im]+'_local_obs2.nc ')
                      os.system(ncap2_str+' -h -O -s  "weights=cos(lat*3.1415/180)" '+cases[im]+'_local_obs2.nc '+cases[im]+'_local_obs3.nc ')
                      os.system(ncap2_str+' -s "SWCF=SWCF*SWCF;" '+cases[im]+'_local_obs3.nc ' +cases[im]+'_local_obs4.nc ')
                      os.system(ncwa_str+' -O  -w weights -a lat,lon  '+cases[im]+'_local_obs4.nc ' +cases[im]+'_local_obs5.nc') 
                      inptrs4 = Dataset(cases[im]+'_local_obs5.nc','r')
                      locals()['MSQ_GLB'] = inptrs4.variables[varisl[ivl]][:]
                      outf.createVariable('MSQ_GLB','f4',('col'))
                      outf.variables['MSQ_GLB'][:]=0
                      outf.variables['MSQ_GLB'][0] = locals()['MSQ_GLB'] #* len(Bnew)

                      inptrs4 = Dataset(str(intll)+'_OBS.nc','r')

                      weigts_sum = 0
                      for slat in range(1, int(edgelat)):
                         for slon in range(1, int(edgelon)):
                             locals()['weights_local'] = inptrs4.variables['weights_'+str(slat)+'_'+str(slon)+'_SWCF'][:]
                             weigts_sum =weigts_sum  + locals()['weights_local']

                      for slat in range(1, int(edgelat)):
                         for slon in range(1, int(edgelon)):
                             tmp_str = ncwa_str+' -O  -w weights -a lat,lon  -d lat,'+str(90.-slat*intll)+','+str(90.+intll-slat*intll)+' -d lon,'+str(1.0*slon*intll-intll)+','+str(1.0*slon*intll)+ ' '+cases[im]+'_local_obs4.nc '+cases[im]+'_local_latlon.nc '
                             os.system(tmp_str)
                             infiles3[im] = cases[im]+'_local_latlon.nc'
                             inptrs3 = Dataset(infiles3[im],'r')
                             locals()['MSQ_'+varisl[ivl]] = inptrs3.variables[varisl[ivl]][:]  #/len(Bnew)
                             locals()['weights_'+str(slat)+'_'+str(slon)] = inptrs4.variables['weights_'+str(slat)+'_'+str(slon)+'_SWCF'][:]
                             outf.createVariable('MSQ_'+varisl[ivl]+'_'+str(slat)+'_'+str(slon),'f4',('col'))
                             outf.variables['MSQ_'+varisl[ivl]+'_'+str(slat)+'_'+str(slon)][:] = 0
                             outf.variables['MSQ_'+varisl[ivl]+'_'+str(slat)+'_'+str(slon)][0] = locals()['MSQ_'+varisl[ivl]]*locals()['weights_'+str(slat)+'_'+str(slon)]/weigts_sum
                             os.system('rm -f '+cases[im]+'_local_latlon.nc ')

                      os.system('rm -f '+cases[im]+'_local_*.nc ')

                 outf.variables[varisl[ivl]+'_GLB'][:]=0
                 outf.variables[varisl[ivl]+'_GLB'][0]=np.sum(locals()[varisl[ivl]][:]*area[:]/np.sum(area))

                 for ir in range(0, len(REG)):
                     locals()[varisl[ivl]+'_'+REG[ir]] = 0
                     locals()['numb_'+REG[ir]] = 0

                 for ir in range(0, len(REG)):
                     for ilat in range(0,nlat):
                         if ((lon[ilat] >= LonW[ir]) & (lon[ilat] < LonE[ir]) & (lat[ilat]>= LatS[ir]) & (lat[ilat] < LatN[ir])):
                            locals()[varisl[ivl]+'_'+REG[ir]] = locals()[varisl[ivl]+'_'+REG[ir]] + locals()[varisl[ivl]][ilat]*area[ilat]
                            locals()['numb_'+REG[ir]] = locals()['numb_'+REG[ir]] + area[ilat]
                     outf.createVariable(varisl[ivl]+'_'+REG[ir],'f4',('col'))
                     outf.variables[varisl[ivl]+'_'+REG[ir]][:] = 0
                     outf.variables[varisl[ivl]+'_'+REG[ir]][0] = locals()[varisl[ivl]+'_'+REG[ir]]/locals()['numb_'+REG[ir]]
                     if (ivl == 0) :
                        outf.createVariable('numb_'+REG[ir],'f4',('col'))
                        outf.variables['numb_'+REG[ir]][:] = 0
                        outf.variables['numb_'+REG[ir]][0] = locals()['numb_'+REG[ir]]/np.sum(area)

                 for slat in range(1, int(edgelat)):
                     for slon in range(1, int(edgelon)):
                         locals()[varisl[ivl]+'_'+str(slat)+'_'+str(slon)] = 0
                         locals()['numb_'+str(slat)+'_'+str(slon)] = 0

                 for slat in range(1, int(edgelat)):
                     for slon in range(1, int(edgelon)):  

                         for ilat in range(0,nlat):
                            if ((lon[ilat] >= 1.0*slon*intll-intll) & (lon[ilat] < 1.0*slon*intll) & (lat[ilat]>= (90.-slat*intll)) & (lat[ilat] < (90.+intll-slat*intll))):
                               locals()[varisl[ivl]+'_'+str(slat)+'_'+str(slon)] = locals()[varisl[ivl]+'_'+str(slat)+'_'+str(slon)] + locals()[varisl[ivl]][ilat]*area[ilat]
                               locals()['numb_'+str(slat)+'_'+str(slon)] = locals()['numb_'+str(slat)+'_'+str(slon)] + area[ilat]
                         outf.createVariable(varisl[ivl]+'_'+str(slat)+'_'+str(slon),'f4',('col'))
                         outf.variables[varisl[ivl]+'_'+str(slat)+'_'+str(slon)][:] = 0
                         outf.variables[varisl[ivl]+'_'+str(slat)+'_'+str(slon)][0] = locals()[varisl[ivl]+'_'+str(slat)+'_'+str(slon)]/locals()['numb_'+str(slat)+'_'+str(slon)]
                         if (ivl == 0) :
                            outf.createVariable('numb_'+str(slat)+'_'+str(slon),'f4',('col'))
                            outf.variables['numb_'+str(slat)+'_'+str(slon)][:] = 0
                            outf.variables['numb_'+str(slat)+'_'+str(slon)][0] = locals()['numb_'+str(slat)+'_'+str(slon)]/np.sum(area)

             for index,parameter in enumerate(list_of_parameters):
                     outf.createVariable(parameter,'f4',('col'))
                     outf.variables[parameter][:]=0
                     outf.variables[parameter][0]=parameter_values[index]
    
             outf.close()

          outf0 = Dataset('./data/'+cseason+'/'+str(intll)+cases[im]+'_Regional.nc','r+')
#          outf0.createVariable(varis[iv]+'_RMSE','f4',('col'))
#          outf0.variables[varis[iv]+'_RMSE'][:]=0
#          outf0.variables[varis[iv]+'_RMSE'][0]=locals()[varis[iv]+'_RMSE_'+cases[im]]
#          outf0.createVariable(varis[iv]+'_RMSEP','f4',('col'))
#          outf0.variables[varis[iv]+'_RMSEP'][:]=0
#          outf0.variables[varis[iv]+'_RMSEP'][0]=locals()[varis[iv]+'_RMSEP_'+cases[im]]
#          outf0.createVariable(varis[iv]+'_RACC','f4',('col'))
#          outf0.variables[varis[iv]+'_RACC'][:]=0
#          outf0.variables[varis[iv]+'_RACC'][0]=locals()[varis[iv]+'_RACC_'+cases[im]]
#          outf0.createVariable(varis[iv]+'_MSE','f4',('col'))
#          outf0.variables[varis[iv]+'_MSE'][:]=0
#          outf0.variables[varis[iv]+'_MSE'][0]=locals()[varis[iv]+'_MSE_'+cases[im]]
          outf0.close()

       grid_x, grid_y = np.mgrid[min(lon):max(lon):100j, min(lat):max(lat):100j]

       A_grid = griddata((lon, lat), A_xy, (grid_x, grid_y), method='cubic')
       grid_shape = (150, 144)
       lon_grid, lat_grid = np.meshgrid(np.linspace(lon.min(), lon.max(), grid_shape[1]),
                                 np.linspace(lat.min(), lat.max(), grid_shape[0]))

       parameter_values=np.zeros(len(list_of_parameters))

       ax = axes[im]
       # grid boxes
       gl = ax.gridlines(
           crs=ccrs.PlateCarree(),
           draw_labels=True,
           linewidth=0.5,
           color='black',
           alpha=0.5,
           linestyle='--'
       )

       #  matplotlib.ticker  FixedLocator 
       gl.xlocator = ticker.FixedLocator(np.arange(-180, 181, intll))
       gl.ylocator = ticker.FixedLocator(np.arange(-90, 91, intll))

# We dont need the labels on top and right
       gl.top_labels   = False
       gl.right_labels = False
       ax.plot = ax.contourf(grid_x, grid_y, A_grid ,
                               levels=levels,  
                               transform=ccrs.PlateCarree(),  # 
                               cmap=cmap1,  # colors
                               extend='both')  

       ax.coastlines()
       ax.set_linewidth=2
       ax.spines['geo'].set_linewidth(1) 
       text_properties = fm.FontProperties( size=14,weight='bold')
       ax.set_title(casenames[im],fontproperties=text_properties,loc='left')


#  END LOP
# OBS 
   ax = axes[im+1]
   gl = ax.gridlines(
       crs=ccrs.PlateCarree(),
       draw_labels=True,
       linewidth=0.5,
       color='black',
       alpha=0.5,
       linestyle='--'
   )

   #  matplotlib.ticker  FixedLocator 
   gl.xlocator = ticker.FixedLocator(np.arange(-180, 181, intll))
   gl.ylocator = ticker.FixedLocator(np.arange(-90, 91, intll))

   gl.top_labels   = False
   gl.right_labels = False
   ax.plot = ax.contourf(lonobs, latobs, B ,
                           levels=levels,
                           transform=ccrs.PlateCarree(),  # 
                           cmap=cmap1,  # colors
                           extend='both')
   ax.coastlines()
   ax.set_title('OBS',fontproperties=text_properties,loc='left')

   if iv == 0:
       lons_mod = [convert_longitude(lons) for lons in lons] 
       ax.scatter( lons_mod, lats, marker='x', s= 7, color='green', label='Markers')
       
       for ir in range(len(REG)):
           LonW_mod = convert_longitude(LonW[ir])
           LonE_mod = convert_longitude(LonE[ir])  #((LonE[ir] - 180) % 360) - 180
           print(LonW_mod)
           xx1, xx2 = LatS[ir], LatN[ir]
           yy1, yy2 = LonW_mod, LonE_mod
           polygon_coords = [(yy1, xx1), (yy2, xx1), (yy2, xx2), (yy1, xx2), (yy1, xx1)]  
           polygon = plt.Polygon(polygon_coords, edgecolor='black', facecolor='white', alpha=0.6)
           ax.add_patch(polygon)
   title_text = f"{varis[iv]}"
#   constrained_layout=True
   cbar_ax = fig.add_axes([0.1, 0.035, 0.8, 0.02])  # [left, bottom, width, height]
   cbar = fig.colorbar(ax.plot, cax=cbar_ax, orientation='horizontal')
   fig.suptitle(title_text,fontsize=16, ha='center', va='center')
   plt.tight_layout()
#   plt.tight_layout(rect=[0, 0, 1, 1])
   plt.savefig(plotname+'.'+ptype, dpi=pixel)
   plt.close()

   del(fig)
   ncdf.close()   
   inptrs.close()
   inptrs2.close()

 return plot2d
