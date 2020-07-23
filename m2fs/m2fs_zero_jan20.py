#!/usr/bin/env python

import numpy as np
import astropy
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
import ccdproc
import astropy.units as u
from astropy.modeling import models
from ccdproc import Combiner
import os
import mycode
import m2fs_process as m2fs
matplotlib.use('TkAgg')

directory='/nfs/nas-0-9/mgwalker.proj/m2fs/'
m2fsrun='jan20' 
datadir=m2fs.get_datadir(m2fsrun)

with open(directory+m2fsrun+'_bias_raw') as f:
    data=f.readlines()[0:]
utdate=[]
file1=[]
file2=[]
for line in data:
    p=line.split()
    utdate.append(str(p[0]))
    file1.append(int(p[1]))
    file2.append(int(p[2]))
utdate=np.array(utdate)
file1=np.array(file1)
file2=np.array(file2)

for ccd in (['b','r']):
    for chip in (['c1','c2','c3','c4']):
        obs_readnoise=[]
        master_processed=[]
        sig_master_processed=[]
        for i in range(0,len(utdate)):
            processed=[]
            sig_processed=[]
            for j in range(file1[i],file2[i]+1):
                filename=datadir+utdate[i]+'/'+ccd+str(j).zfill(4)+chip+'.fits'
                
                data=astropy.nddata.CCDData.read(filename,unit=u.adu)#header is in data.meta
                print(filename,data.header['object'],data.header['binning'])

                oscan_subtracted=ccdproc.subtract_overscan(data,overscan=data[:,1024:],overscan_axis=1,model=models.Polynomial1D(3),add_keyword={'oscan_corr':'Done'})
                trimmed1=ccdproc.trim_image(oscan_subtracted[:,:1024],add_keyword={'trim1':'Done'})
                trimmed2=ccdproc.trim_image(trimmed1[:1028,:1024],add_keyword={'trim2':'Done'})
                array1d=trimmed2.data.flatten()
                gain=np.float(trimmed2.header['egain'])
                keep=np.where(np.abs(array1d)<100.)[0]#remove crazy outliers
                obs_readnoise.append(np.std(array1d[keep]*gain))
#                data_with_deviation=ccdproc.create_deviation(trimmed2,gain=data.meta['egain']*u.electron/u.adu,readnoise=data.meta['enoise']*u.electron)
#                gain_corrected=ccdproc.gain_correct(data_with_deviation,data_with_deviation.meta['egain']*u.electron/u.adu,add_keyword={'gain_corr':'Done'})
#                cr_cleaned=ccdproc.cosmicray_lacosmic(trimmed2,sigclip=5,gain_apply=False,gain=0.68,readnoise=2.7)
#                sig_cr_cleaned=cr_cleaned.uncertainty._array

#                data_with_deviation=ccdproc.create_deviation(data,gain=data.meta['egain']*u.electron/u.adu,readnoise=data.meta['enoise']*u.electron)
#                gain_corrected=ccdproc.gain_correct(data_with_deviation,data_with_deviation.meta['egain']*u.electron/u.adu,add_keyword={'gain_corr':'Done'})
#                cr_cleaned=ccdproc.cosmicray_lacosmic(gain_corrected,sigclip=5)
#                oscan_subtracted=ccdproc.subtract_overscan(cr_cleaned,overscan=cr_cleaned[:,1024:],overscan_axis=1,model=models.Polynomial1D(3),add_keyword={'oscan_corr':'Done'})
#                trimmed1=ccdproc.trim_image(oscan_subtracted[:,:1024],add_keyword={'trim1':'Done'})
#                trimmed2=ccdproc.trim_image(trimmed1[:1028,:1024],add_keyword={'trim2':'Done'})
#                trimmed.append(trimmed2)
#                sig_trimmed.append(trimmed2.uncertainty._array)
                processed.append(trimmed2)
#                sig_processed.append(sig_cr_cleaned)
                master_processed.append(trimmed2)
#                sig_master_processed.append(sig_cr_cleaned)
#            processed=np.array(processed)
#            sig_processed=np.array(sig_processed)

#            c=Combiner(processed)
#            c.clip_extrema(nlow=1,nhigh=1)
#            old_n_masked=0
#            new_n_masked=c.data_arr.mask.sum()
#            while (new_n_masked > old_n_masked):
#                c.sigma_clipping(low_thresh=3,high_thresh=3,func=np.ma.median)
#                old_n_masked=new_n_masked
#                new_n_masked=c.data_arr.mask.sum()
#            ccdall=c.average_combine()
#            ccdall.write(directory+m2fsrun+'_'+ccd+'_'+chip+'_master_bias'+str(i+1)+'.fits',overwrite=True)
        obs_readnoise=np.array(obs_readnoise)
        c=Combiner(master_processed)
        c.clip_extrema(nlow=1,nhigh=1)
        old_n_masked=0
        new_n_masked=c.data_arr.mask.sum()
        while (new_n_masked > old_n_masked):
#            c.sigma_clipping(func=np.ma.median)
            c.sigma_clipping(low_thresh=3,high_thresh=3,func=np.ma.median)
            old_n_masked=new_n_masked
            new_n_masked=c.data_arr.mask.sum()
#        c.clip_extrema(nlow=5,nhigh=5)
#        c.weights=1./sig_master_bias**2
#        ccdall=c.average_combine(uncertainty_func=mycode.stdmean)
        ccdall=c.average_combine()
        ccdall[0].header['obs_rdnoise']=str(np.median(obs_readnoise))
        ccdall.write(directory+m2fsrun+'_'+ccd+'_'+chip+'_master_bias.fits',overwrite=True)

