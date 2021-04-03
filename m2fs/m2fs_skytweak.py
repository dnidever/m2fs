#!/usr/bin/env python

import os
import numpy as np
import dill as pickle
from dlnpyutils import bindata, utils as dln
from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('MacOSX')

# take flux, subtract scaled sky spectrum, then fit a continuum and subtract it
# then check the median (subtracted) flux of the pixels with bright sky lines
# want to minimize that

import astropy.units as u
from astropy.modeling import models,fitting
from copy import deepcopy
import copy


#continuum_init = models.Chebyshev1D(degree=continuum_rejection_order)
#fitter = fitting.LinearLSQFitter()
#lamb = spec1d.spectral_axis
#y = np.ma.masked_array(spec1d.flux,mask=spec1d.mask)
#continuum=fitter(continuum_init,lamb.value[y.mask==False],y[y.mask==False])

def func_sky(x,pars):
    pass

    
def fitsky(ext,sky):

    out = curve_fit('func_sky',ext.spec1d_flux.data,ext.spec1d_uncertainty.quantity.val)

def get_continuum(flux):
    # Fit continuum
    x = np.arange(len(flux))
    ybin,bin_edges,binnumber = bindata.binned_statistic(x,flux,statistic='median',bins=10)
    xbin = bin_edges[0:-1]+(bin_edges[1]-bin_edges[0])*0.5
    cont = dln.interp(xbin,ybin,x,extrapolate=True)
    return cont
    
def skyresid(ext,sky,scl,mask):
    resid = ext-sky*scl
    # Fit continuum and subtract
    cont = get_continuum(resid)
    resid -= cont
    return np.median(resid[mask])
        
def skytweak(extfile,skyfile,plugfile,diag=False):
    """ Tweak the Sky subtraction."""

    #extfile = 'ut20131123/b0236_extract1d_array.pickle'
    #skyfile = 'ut20131123/b0236_sky_array.pickle'
    #plugfile = 'ut20131123/b0236_plugmap.pickle'

    print('Tweaking sky subtraction')
    print('extract1d file = '+extfile)
    print('sky file = '+skyfile)
    print('plugmap file = '+plugfile)    
    
    if os.path.exists(extfile)==False:
        raise ValueError(extfile+' NOT FOUND')
    if os.path.exists(skyfile)==False:
        raise ValueError(skyfile+' NOT FOUND')
    if os.path.exists(plugfile)==False:
        raise ValueError(plugfile+' NOT FOUND')    
    
    ext = pickle.load(open(extfile,'rb'))    
    sky = pickle.load(open(skyfile,'rb'))
    plugmap = pickle.load(open(plugfile,'rb'))

    #skysub = pickle.load(open('ut20131123/b0236_skysubtract_array.pickle','rb'))
    #ext1d = pickle.load(open('ut20131123/b0236_extract1d_array.pickle','rb'))
    #ext = pickle.load(open('ut20131123/b0236_throughputcorr_array.pickle','rb'))
    #thru = pickle.load(open('ut20131123/b0236_throughput_array.pickle','rb'))
    #plugmap = pickle.load(open('ut20131123/b0236_plugmap.pickle','rb')) 
    
    nap = len(ext)
    outsky = []
    outskysub = []
    outscale = []
    # Loop over apertures
    for i in range(nap):
        ext1 = ext[i]
        sky1 = sky[i]
        if np.nansum(sky1.spec1d_flux.value)>0:
            skycont = get_continuum(sky1.spec1d_flux.value)
            mask = ((sky1.spec1d_flux.value-skycont) > 25) & np.isfinite(ext1.spec1d_flux.value) & np.isfinite(sky1.spec1d_flux.value)
        
            scaling = np.linspace(0.2,2.0,50)
            resid = np.zeros(len(scaling),float)
            for j in range(len(scaling)):
                resid[j] = skyresid(ext1.spec1d_flux.value,sky1.spec1d_flux.value,scaling[j],mask)
            scaling2 = np.linspace(0.2,2.0,1000)
            resid2 = dln.interp(scaling,resid,scaling2)
            bestind = np.argmin(np.abs(resid2))
            bestscale = scaling2[bestind]
            print(str(i+1)+' '+str(bestscale))
            outscale.append(bestscale)
            outsky1 = deepcopy(sky1)
            outsky1.spec1d_flux *= bestscale
            outsky1.spec1d_uncertainty.array *= bestscale            
            outsky.append(outsky1)
            outskysub1 = deepcopy(ext1)
            outskysub1.spec1d_flux -= outsky1.spec1d_flux
            outskysub.append(outskysub1)

            eflux = ext1.spec1d_flux.value.copy()
            sflux = sky1.spec1d_flux.value.copy()
            bestresid = eflux-sflux*bestscale
            # Fit continuum and subtract
            bestcont = get_continuum(bestresid)
    
            if diag==True:
                eflux[(~np.isfinite(eflux)) | (ext1.spec1d_uncertainty.quantity.value>900)] = np.nan                
                fig = plt.figure(figsize=(10,8))
                #fig,ax = plt.subplots()
                #fig.set_figheight(8)
                #fig.set_figwidth(10)
                plt.plot(eflux)
                plt.plot(bestcont)
                plt.plot(outsky1.spec1d_flux.value+bestcont,alpha=0.7)
                plt.plot([0,2048],[0,0],linestyle='--',c='gray',alpha=0.8)
                #plt.plot(outsky1.spec1d_flux)
                plt.plot(outskysub1.spec1d_flux.value-150)
                plt.plot(bestcont-150)
                plt.plot([0,2048],[-150,-150],linestyle='--',c='gray',alpha=0.8)
                plt.xlabel('X')
                plt.ylabel('Flux')
                plt.xlim(0,2048)
                plt.ylim(-400,800)
                plt.title('Aperture '+str(ext1.aperture)+' Scale=%5.2f' % bestscale)
                plt.show()
        else:
            print(str(i+1)+' skipping')
            outscale.append(1.0)
            outsky.append(sky1)
            outskysub.append(ext1)
            
    return outscale,outsky,outskysub


if __name__ == "__main__":
    parser = ArgumentParser(description='Run Doppler fitting on spectra')
    parser.add_argument('extfile', type=str, nargs=1, help='Extract1D filename')
    parser.add_argument('--skyfile', type=str, nargs=1, default='', help='Skyfile')
    parser.add_argument('--plugfile', type=str, nargs=1, default='', help='Plugmap file')
    parser.add_argument('--outfile', type=str, nargs=1, default='', help='Output filename')
    parser.add_argument('-d','--diag', action='store_true', help='Show diagnostic information')
    args = parser.parse_args()
    extfile = args.extfile[0]
    skyfile = dln.first_el(args.skyfile)
    plugfile = dln.first_el(args.plugfile)    
    outfile = dln.first_el(args.outfile)
    diag = args.diag

    datadir = os.path.dirname(extfile)
    base = os.path.basename(extfile)
    base = base.split('_')[0]
    if skyfile=='':
        skyfile = datadir+'/'+base+'_sky_array.pickle'
    if plugfile=='':
        plugfile = datadir+'/'+base+'_plugmap.pickle'
    if outfile=='':
        outfile = datadir+'/'+base+'_skytweak.pickle'        
    
        
    #extfile = 'ut20131123/b0236_extract1d_array.pickle'
    #skyfile = 'ut20131123/b0236_sky_array.pickle'
    #plugfile = 'ut20131123/b0236_plugmap.pickle'
    #skytweakfile = 'skytweak.pickle'
    
    scale,sky,skysub = skytweak(extfile,skyfile,plugfile,diag=diag)

    print('Writing results to '+outfile)
    pickle.dump((scale,sky,skysub),open(outfile,'wb'))
    # Write new sky_array.pickle and skysubtract_array.pickle files
    outskyfile = datadir+'/'+base+'_sky_tweak_array.pickle'
    print('Writing new sky to '+outskyfile)
    if os.path.exists(outskyfile): os.remove(outskyfile)
    pickle.dump(sky,open(outskyfile,'wb'))
    outskysubfile = datadir+'/'+base+'_skysubtract_tweak_array.pickle'
    print('Writing new skysub to '+outskysubfile)
    if os.path.exists(outskysubfile): os.remove(outskysubfile)
    pickle.dump(skysub,open(outskysubfile,'wb'))    
