#!/usr/bin/env python

import astropy
import dill as pickle
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
from astropy.nddata import Cutout2D
from astropy.nddata import NDDataRef
import astropy.units as u
from astropy.modeling import models
import m2fs_process as m2fs
import os
from os import path
import utils
import specutils
from specutils.spectra import Spectrum1D
import scipy
from astropy.nddata import StdDevUncertainty
from astropy.visualization import quantity_support
import numpy as np
from specutils.fitting import fit_lines
from specutils.analysis import centroid
from specutils.analysis import fwhm
from specutils.analysis import line_flux
from specutils.analysis import equivalent_width
from specutils import SpectralRegion
from specutils.manipulation import extract_region
from specutils.fitting import fit_generic_continuum
from astropy.modeling import models,fitting
from specutils.fitting import estimate_line_parameters
from ccdproc import Combiner
from scipy import interpolate
from astropy import time, coordinates as coord, units as u
#from astropy.coordinates import SkyCoord, EarthLocation
from mpl_toolkits.axes_grid1 import make_axes_locatable
#matplotlib.use('pdf')
matplotlib.use('TkAgg')

directory = '/nfs/nas-0-9/mgwalker.proj/m2fs/'

with open(directory+'m2fs_zero0.py') as f:
    zero_text = f.readlines()[0:]
with open(directory+'m2fs_dark0.py') as f:
    dark_text = f.readlines()[0:]
with open(directory+'m2fs_reduce0.py') as f:
    reduce_text = f.readlines()[0:]
with open(directory+'m2fs_apertures0.py') as f:
    apertures_text = f.readlines()[0:]
with open(directory+'m2fs_apertures_noflat0.py') as f:
    apertures_noflat_text = f.readlines()[0:]

m2fsrun0 = ['feb14','dec14','feb15','jul15','sep15','nov15','feb16','jun16','aug16','nov16','feb17','may17',
            'sep17','nov17','feb18','may18','aug18','nov18','feb19','may19','aug19','nov19','jan20']
#m2fsrun0 = ['aug16']

for i in range(0,len(m2fsrun0)):
    zero_out = 'm2fs_zero_'+m2fsrun0[i]+'.py'
    dark_out = 'm2fs_dark.py'
    reduce_out = 'm2fs_reduce_'+m2fsrun0[i]+'.py'
    apertures_out = 'm2fs_apertures_'+m2fsrun0[i]+'.py'
    apertures_noflat_out = 'm2fs_apertures_noflat_'+m2fsrun0[i]+'.py'
    zero_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_zero_'+m2fsrun0[i]+'.slurm'
    dark_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_dark.slurm'#only one master dark for all runs
    reduce_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_reduce_'+m2fsrun0[i]+'.slurm'
    apertures_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_apertures_'+m2fsrun0[i]+'.slurm'
    apertures_slurm_noflat_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_apertures_noflat_'+m2fsrun0[i]+'.slurm'

    g1 = open(zero_out,'w')
    for line in zero_text:
        if 'enter_run_here' in line:
            g1.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        else:
            g1.write(line)
    g1.close()

    g1 = open(dark_out,'w')
    for line in dark_text:
#        if 'enter_run_here' in line:
#            g1.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
#        else:
        g1.write(line)
    g1.close()

    g1 = open(reduce_out,'w')
    for line in reduce_text:
        if 'enter_run_here' in line:
            g1.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        else:
            g1.write(line)
    g1.close()

    g1 = open(apertures_out,'w')
    for line in apertures_text:
        if 'enter_run_here' in line:
            g1.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        else:
            g1.write(line)
    g1.close()

    g1 = open(apertures_noflat_out,'w')
    for line in apertures_noflat_text:
        if 'enter_run_here' in line:
            g1.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        else:
            g1.write(line)
    g1.close()

    g1 = open(zero_slurm_out,'w')
    g1.write('#!/bin/bash \n')
#    g1.write('#SBATCH -N 1 \n')
    g1.write('#SBATCH --ntasks=1 \n')
    g1.write('#SBATCH --time=0-07:59 \n')
#    g1.write('#SBATCH --partition=short \n')
    g1.write('#SBATCH --mem=100000 \n')
    g1.write('#SBATCH -o m2fs_zero_'+m2fsrun0[i]+'.o \n')
    g1.write('#SBATCH -e m2fs_zero_'+m2fsrun0[i]+'.err \n')
    g1.write(' \n')
    g1.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+zero_out+' \n')
    g1.close()

    g1 = open(dark_slurm_out,'w')
    g1.write('#!/bin/bash \n')
#    g1.write('#SBATCH -N 1 \n')
    g1.write('#SBATCH --ntasks=1 \n')
    g1.write('#SBATCH --time=0-07:59 \n')
#    g1.write('#SBATCH --partition=short \n')
    g1.write('#SBATCH --mem=100000 \n')
    g1.write('#SBATCH -o m2fs_dark.o \n')
    g1.write('#SBATCH -e m2fs_dark.err \n')
    g1.write(' \n')
    g1.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+dark_out+' \n')
    g1.close()

    g1 = open(reduce_slurm_out,'w')
    g1.write('#!/bin/bash \n')
#    g1.write('#SBATCH -N 1 \n')
    g1.write('#SBATCH --ntasks=1 \n')
    g1.write('#SBATCH --time 0-07:59 \n')
#    g1.write('#SBATCH --partition=short \n')
    g1.write('#SBATCH --mem=100000 \n')
    g1.write('#SBATCH -o m2fs_reduce_'+m2fsrun0[i]+'.o \n')
    g1.write('#SBATCH -e m2fs_reduce_'+m2fsrun0[i]+'.err \n')
    g1.write(' \n')
    g1.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+reduce_out+' \n')
    g1.close()

    g1 = open(apertures_slurm_out,'w')
    g1.write('#!/bin/bash \n')
#    g1.write('#SBATCH -N 1 \n')
    g1.write('#SBATCH --ntasks=1 \n')
    g1.write('#SBATCH --time=0-07:59 \n')
#    g1.write('#SBATCH --partition=short \n')
    g1.write('#SBATCH --mem=100000 \n')
    g1.write('#SBATCH -o m2fs_apertures_'+m2fsrun0[i]+'.o \n')
    g1.write('#SBATCH -e m2fs_apertures_'+m2fsrun0[i]+'.err \n')
    g1.write(' \n')
    g1.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+apertures_out+' \n')
    g1.close()

    g1 = open(apertures_slurm_noflat_out,'w')
    g1.write('#!/bin/bash \n')
#    g1.write('#SBATCH -N 1 \n')
    g1.write('#SBATCH --ntasks=1 \n')
    g1.write('#SBATCH --time=0-07:59 \n')
    g1.write('#SBATCH --partition=long \n')
    g1.write('#SBATCH --mem=100000 \n')
    g1.write('#SBATCH -o m2fs_apertures_noflat_'+m2fsrun0[i]+'.o \n')
    g1.write('#SBATCH -e m2fs_apertures_noflat_'+m2fsrun0[i]+'.err \n')
    g1.write(' \n')
    g1.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+apertures_noflat_out+' \n')
    g1.close()

    utdate = []
    file1 = []
    file2 = []
    flatfile = []
    tharfile = []
    field_name = []
    scifile = []
    fibermap_file = []
    fiber_changes = []
    obj = []

    with open(directory+m2fsrun0[i]+'_science_raw') as f:
        data = f.readlines()[0:]
    for line in data:
        p = line.split()
        if p[0]!='none':
            utdate.append(str(p[0]))
            file1.append(int(p[1]))
            file2.append(int(p[2]))
            flatfile.append(p[3])
            tharfile.append(p[4])
            field_name.append(p[5])
            scifile.append(p[6])
            fibermap_file.append(p[7])
            fiber_changes.append(p[8])
            obj.append(p[9])
    utdate = np.array(utdate)
    file1 = np.array(file1)
    file2 = np.array(file2)
    flatfile = np.array(flatfile)
    tharfile = np.array(tharfile)
    field_name = np.array(field_name)
    scifile = np.array(scifile)
    fibermap_file = np.array(fibermap_file)
    fiber_changes = np.array(fiber_changes)
    obj = np.array(obj)

    for j in range(0,len(utdate)):
        science_raw_out = m2fsrun0[i]+'_'+str(j+1)+'_science_raw'
        apertures2_out = 'm2fs_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py'
        apertures2_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.slurm'
        apertures_noflat2_out = 'm2fs_apertures_noflat_'+m2fsrun0[i]+'_'+str(j+1)+'.py'
        apertures_noflat2_slurm_out = '/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_apertures_noflat_'+m2fsrun0[i]+'_'+str(j+1)+'.slurm'
        g1 = open(science_raw_out,'w')
        g2 = open(apertures2_out,'w')
        g3 = open(apertures2_slurm_out,'w')
        g4 = open(apertures_noflat2_out,'w')
        g5 = open(apertures_noflat2_slurm_out,'w')
        g1.write(utdate[j]+' '+str(file1[j]).zfill(4)+' '+str(file2[j]).zfill(4)+' '+flatfile[j]+' '+tharfile[j]+' '+field_name[j]+' '+scifile[j]+' '+fibermap_file[j]+' '+fiber_changes[j]+' '+obj[j]+' \n')
        g1.close()
        
        for line in apertures_text:
            if 'enter_run_here' in line:
                g2.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
            elif 'with open(directory+m2fsrun+\'_science_raw\') as f:' in line:
                g2.write('with open(directory+m2fsrun+\''+'_'+str(j+1)+'_science_raw\') as f: \n')
            else:
                g2.write(line)
        g2.close()

        g3.write('#!/bin/bash \n')
#        g3.write('#SBATCH -N 1 \n')
        g3.write('#SBATCH --ntasks=1 \n')
        g3.write('#SBATCH --time=0-07:59 \n')
#        g3.write('#SBATCH --partition=short \n')
        g3.write('#SBATCH --mem=100000 \n')
        g3.write('#SBATCH -o m2fs_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.o \n')
        g3.write('#SBATCH -e m2fs_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.err \n')
        g3.write(' \n')
        g3.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/m2fs_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py \n')
        g3.close()

        for line in apertures_noflat_text:
            if 'enter_run_here' in line:
                g4.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
            elif 'with open(directory+m2fsrun+\'_science_raw\') as f:' in line:
                g4.write('with open(directory+m2fsrun+\''+'_'+str(j+1)+'_science_raw\') as f: \n')
            else:
                g4.write(line)
        g4.close()

        g5.write('#!/bin/bash \n')
#        g5.write('#SBATCH -N 1 \n')
        g5.write('#SBATCH --ntasks=1 \n')
        g5.write('#SBATCH --time=0-07:59 \n')
        g5.write('#SBATCH --partition=long \n')
        g5.write('#SBATCH --mem=100000 \n')
        g5.write('#SBATCH -o m2fs_apertures_noflat_'+m2fsrun0[i]+'_'+str(j+1)+'.o \n')
        g5.write('#SBATCH -e m2fs_apertures_noflat_'+m2fsrun0[i]+'_'+str(j+1)+'.err \n')
        g5.write(' \n')
        g5.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/m2fs_apertures_noflat_'+m2fsrun0[i]+'_'+str(j+1)+'.py \n')
        g5.close()
