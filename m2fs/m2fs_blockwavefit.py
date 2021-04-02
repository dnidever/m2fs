#!/usr/bin/env python

import os
import numpy as np
import dill as pickle
from astropy.io import fits
from astropy.modeling import models,fitting
from dlnpyutils import bindata, utils as dln
from argparse import ArgumentParser
from doppler.rv import ccorrelate
from scipy.interpolate import BSpline
from scipy.interpolate import UnivariateSpline,LSQUnivariateSpline

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
import m2fs_process as m2fs

def getaplines(idlines,apinfo,ap,diag=False):

    #aperture = [id.aperture for id in idlines]
    #aperture = np.array(aperture)
    apind, = np.where(apinfo['aperture']==ap)
    apind = apind[0]
    idline1 = idlines[apind]

    goodap, = np.where(apinfo['good']==True)
    
    # This is a "good" aperture
    if apind in goodap:
        #print(str(ap)+' is a good aperture')
        nlines = len(idline1.fit_lines.fit)
        linecat = np.zeros(nlines,dtype=np.dtype([('aperture',int),('num',int),('height',float),('center',float),
                                                  ('sigma',float),('wave0',float),('wave',float),('match',int),('outlier',bool)]))
        linecat['aperture'][:] = idline1.aperture
        linecat['num'][:] = np.arange(nlines)+1
        for j in range(nlines):
            linecat['height'][j] = idline1.fit_lines.fit[j].parameters[0]
            linecat['center'][j] = idline1.fit_lines.fit[j].parameters[1]
            linecat['sigma'][j] = idline1.fit_lines.fit[j].parameters[2]
            wav = idline1.wav[j]
            if isinstance(wav,float):
                linecat['wave'][j] = wav
                linecat['match'][j] = 1
        return linecat

    print(str(ap)+' is a bad aperture')
    
    # This is a "bad" aperture
    spec = extract1d_array[apind]
    # Find a neighboring good aperture
    bestind = goodap[np.argmin(np.abs(apinfo['aperture'][apind]-apinfo['aperture'][goodap]))]
    #if idline1.aperture==61:
    #    bestind = 62
    print('Best reference aperture: '+str(apinfo['aperture'][bestind]))
    refap = idlines[bestind]
    refspec = extract1d_array[bestind]
    # Cross-correlate the spectra to get shift
    maxlag = 10
    lag = np.arange(2*maxlag+1)-maxlag
    rflux = refspec.spec1d_flux.value
    rflux[refspec.spec1d_uncertainty.quantity.value>900] = np.nan
    sflux = spec.spec1d_flux.value
    sflux[spec.spec1d_uncertainty.quantity.value>900] = np.nan        
    # cross-correlate chunks
    nchunks = 5
    xshift = np.zeros(nchunks,float)
    mnx = np.zeros(nchunks,float)
    nbin = npix//nchunks
    for k in range(nchunks):
        mnx[k] = np.mean([k*nbin,(k+1)*nbin-1])
        cc = ccorrelate(rflux[k*nbin:(k+1)*nbin],sflux[k*nbin:(k+1)*nbin],lag)
        lag2 = np.linspace(-maxlag,maxlag,1000)
        cc2 = dln.interp(lag,cc,lag2)
        bestshift = lag2[np.argmax(cc2)]
        xshift[k] = bestshift
    print('shifts = '+str(xshift))
    cc = ccorrelate(rflux,sflux,lag)
    lag2 = np.linspace(-maxlag,maxlag,1000)
    cc2 = dln.interp(lag,cc,lag2)
    bestshift = lag2[np.argmax(cc2)]        
    print('best shift = '+str(bestshift))

    #if ap==61:
    #    #plt.plot(x,rflux)
    #    #plt.plot(x+bestshift,sflux)
    #    #plt.show()
    #    #import pdb; pdb.set_trace()
    #    print('KLUDGE for aperture 61')
    #    bestshift = 0.0
        
    if diag==True:
        plt.plot(x,rflux)
        plt.plot(x-bestshift,sflux)
        plt.show()
    
    # Shift reference wavelength array
    refwav = refap.func(x)
    drefwav = refwav[1:]-refwav[0:-1]
    drefwav = np.append(drefwav,drefwav[-1])
    xshiftall = dln.interp(mnx,xshift,np.arange(npix))
    #wav0 = refwav.copy() - drefwav*xshiftall
    wav0 = refwav.copy() - drefwav*bestshift        
    # Get initial wavelengths for the lines
    nlines = len(idline1.fit_lines.fit)
    linecat = np.zeros(nlines,dtype=np.dtype([('aperture',int),('num',int),('height',float),('center',float),
                                              ('sigma',float),('wave0',float),('wave',float),('match',int),('outlier',bool)]))
    linecat['aperture'][:] = idline1.aperture
    linecat['num'][:] = np.arange(nlines)+1
    for j in range(nlines):
        linecat['height'][j] = idline1.fit_lines.fit[j].parameters[0]
        linecat['center'][j] = idline1.fit_lines.fit[j].parameters[1]
        linecat['sigma'][j] = idline1.fit_lines.fit[j].parameters[2]            
    linecat['wave0'][:] = dln.interp(x,wav0,linecat['center'])
    # Try to ID them using the template lines
    match = np.zeros(ntlines,int)-1
    dwave = np.zeros(ntlines,float)-1
    for j in range(ntlines):
        bestline = np.argmin(np.abs(tlinecat['wave'][j]-linecat['wave0']))
        maxdiff = 0.5  #1.0
        if np.abs(tlinecat['wave'][j]-linecat['wave0'][bestline])<maxdiff:
            match[j] = bestline
            dwave[j] = tlinecat['wave'][j]-linecat['wave0'][bestline]
    gmatch, = np.where(match>-1)
    linecat['wave'][match[gmatch]] = tlinecat['wave'][gmatch]
    linecat['match'][match[gmatch]] = 1

    return linecat

if __name__ == "__main__":
    parser = ArgumentParser(description='Run Doppler fitting on spectra')
    parser.add_argument('idlinefile', type=str, nargs=1, help='M2FS ID lines translate filename')
    parser.add_argument('--outfile', type=str, nargs=1, default='', help='Output file')
    parser.add_argument('--directory', type=str, nargs=1, default='./', help='Main directory with the data')
    parser.add_argument('-d','--diag', action='store_true', help='Show diagnostic information')
    args = parser.parse_args()
    idlinefile = args.idlinefile[0]
    outfile = dln.first_el(args.outfile)
    directory = dln.first_el(args.directory)
    diag = args.diag
    
    npix = 2048
    x = np.arange(npix)
    
    print('Fitting block wavelength solutions for '+idlinefile)

    # Get template name
    with open(directory+'arc_templates') as f:
        data = f.readlines()
    arcfilename = []
    arcfiltername = []
    for line in data:
        p = line.split()
        arcfilename.append(p[0])
        arcfiltername.append(p[1])
    arcfilename = np.array(arcfilename)
    arcfiltername = np.array(arcfiltername)

    datadir = os.path.dirname(idlinefile)
    root0 = os.path.basename(idlinefile)
    root0 = root0.split('_')[0]
    
    header = fits.getheader(datadir+'/'+root0+'_stitched.fits')
    filtername=header['FILTER']
    if filtername=='Mgb_Rev2':
        filtername='Mgb_HiRes'
    if filtername=='CalRT_O41' and header['FF-THNE'] > 0:
        filtername='CaT_ThNe'
    if filtername=='CalRT_O41' and header['FF-THAR'] > 0:
        filtername='CaT_ThAr'
    if filtername=='CalRT_O41' and header['FF-THAR']==0 and header['EXPTYPE']=='Object':
        filtername='CaT_Sky'

    # Load the template information
    arc = np.where(arcfiltername==filtername)[0][0]
    template_root = str(arcfilename[arc])
    extract1d_array_template_file = template_root+'_extract1d_array.pickle'
    id_lines_template_file = template_root+'_id_lines_template.pickle'
    extract1d_array_template = pickle.load(open(extract1d_array_template_file,'rb'))
    id_lines_template = pickle.load(open(id_lines_template_file,'rb'))
    
    ntlines = len(id_lines_template.fit_lines.fit)
    tlinecat = np.zeros(ntlines,dtype=np.dtype([('num',int),('height',float),('center',float),('sigma',float),('wave',float)]))
    tlinecat['num'][:] = np.arange(ntlines)+1
    for j in range(ntlines):
            tlinecat['height'][j] = id_lines_template.fit_lines.fit[j].parameters[0]
            tlinecat['center'][j] = id_lines_template.fit_lines.fit[j].parameters[1]
            tlinecat['sigma'][j] = id_lines_template.fit_lines.fit[j].parameters[2]
            wav = id_lines_template.wav[j]
            if isinstance(wav,float):
                tlinecat['wave'][j] = wav
    g, = np.where(tlinecat['wave'] > 10)  # only keep lines with wavelengths
    tlinecat = tlinecat[g]
    ntlines = len(tlinecat)
    
    # Load the arc data
    idlines = pickle.load(open(idlinefile,'rb'))
    n = len(idlines)
    extract1d_array_file = datadir+'/'+root0+'_extract1d_array.pickle'
    extract1d_array = pickle.load(open(extract1d_array_file,'rb'))
    
    # Gather RMS and Npoints
    apinfo = np.zeros(n,dtype=np.dtype([('aperture',int),('npoints',int),('rms',float),('good',bool)]))
    for i in range(n):
        apinfo['aperture'][i] = idlines[i].aperture
        apinfo['rms'][i] = idlines[i].rms
        apinfo['npoints'][i] = idlines[i].npoints                
    used, = np.where(apinfo['npoints']>0)
    medrms = np.median(apinfo['rms'][used])
    mednpoints = np.median(apinfo['npoints'][used])

    # Bad ones to fix
    #badap1, = np.where((npoints>0) & ((npoints<0.5*mednpoints) | (rms>3*medrms)))
    goodap1, = np.where((apinfo['npoints']>0.5*mednpoints) & (apinfo['rms']<1.5*medrms))
    
    # Compare wavelenths to median wavelength plus shift
    goodwav0 = np.zeros((npix,len(goodap1)),float)
    for i in range(len(goodap1)):
        goodwav0[:,i] = idlines[goodap1[i]].func(x)
    medwav0 = np.median(goodwav0,axis=1)
    goodwav = np.zeros((npix,len(goodap1)),float)    
    for i in range(len(goodap1)):
        wav1 = idlines[goodap1[i]].func(x)
        diff = wav1-medwav0
        meddiff = np.median(diff)
        goodwav[:,i] = wav1-meddiff
    medwav = np.median(goodwav,axis=1)

    # Second cut at good and bad apertures by comparing to median wave array
    sig = np.zeros(128,float)+999999.0
    std = np.zeros(128,float)+999999.0
    for i in range(128):
        if idlines[i].npoints>0:
            wav1 = idlines[i].func(x)
            std[i] = np.std(wav1-medwav)
            sig[i] = dln.mad(wav1-medwav)            


    #goodap, = np.where((npoints>0) & (sig<1.0) & (std<1.0) & (npoints>0.5*mednpoints))
    #badap, = np.where((npoints>0) & ((sig>=1.0) | (std>=1.0) | (npoints<=0.5*mednpoints)))
    goodap, = np.where((apinfo['npoints']>0) & (sig<1.0) & (std<1.0) & (apinfo['rms']<0.1) & (apinfo['npoints']>20))
    badap, = np.where((apinfo['npoints']>0) & ((sig>=1.0) | (std>=1.0) | (apinfo['rms']>=0.1) | (apinfo['npoints']<=20)))
    nbad = len(badap)
    print(str(nbad)+' bad apertures to fix')
    apinfo['good'][goodap] = True

    
    #goodap, = np.where((npoints>0) & (sig<1.0) & (std<1.0) & (npoints>0.5*mednpoints) & (rms<0.08))
    #badap, = np.where(rms > 0.08)
    #nbad = len(badap)
    #import pdb; pdb.set_trace()

    # Loop over blocks
    nblocks = 8
    order = 5 # 4
    for b in range(nblocks):
        print('Fitting wavelength solution for block '+str(b+1))
        # apertures with some lines
        gdapind, = np.where((apinfo['npoints']>0) & (apinfo['aperture']>=(b*16+1)) & (apinfo['aperture']<((b+1)*16+1)))
        gdap = apinfo['aperture'][gdapind]
        nap = len(gdap)
        print(str(len(gdap))+' apertures')
        
        # Get aperture lines
        linecat = None
        for ap in gdap:
            linecat1 = getaplines(idlines,apinfo,ap,diag=diag)
            print(ap,len(linecat1))
            if linecat is None:
                linecat = linecat1
            else:
                linecat = np.hstack((linecat,linecat1))
        nlinecat = len(linecat)
        
        # Fit 4th or 5th order polynomial plus a separate offset for each aperture

        # 1) Initial polynomial fit to all lines in all fibers
        gline, = np.where(linecat['match']==1)
        func_init = models.Legendre1D(degree=order,domain=[0,npix])
        fitter = fitting.LinearLSQFitter()
        func1 = fitter(func_init,linecat['center'][gline],linecat['wave'][gline])
        dwave1 = func1(linecat['center'][gline]) - linecat['wave'][gline]
        #coef1 = dln.poly_fit(linecat['center'][gline],linecat['wave'][gline],4)
        #dwave1 = dln.poly(linecat['center'][gline],coef1) - linecat['wave'][gline]
        sigwave1 = dln.mad(dwave1)
        print('Initial fit coefficients ',func1.parameters)
        #print('Initial fit coefficients ',coef1)        
        print('Initial sigma = %f6.2' % sigwave1)

        #xx = linecat['center'][gline]
        #yy = linecat['wave'][gline]
        #si = np.argsort(xx)
        #spl = UnivariateSpline(xx[si],yy[si],k=3,s=10)
        #dwave1 = spl(xx) - linecat['wave'][gline]

        #import matplotlib.pyplot as plt
        #import matplotlib
        #matplotlib.use('MacOSX')
        #plt.scatter(xx,dwave1,c=linecat['aperture'][gline])
        #plt.show()
        
        #import pdb; pdb.set_trace()
        
        
        # 2) Initial offsets
        woff1 = np.zeros(nap,float)
        wavecorr1 = np.zeros(nlinecat,float)
        for i,ap in enumerate(gdap):
            ind, = np.where((linecat['aperture']==ap) & (linecat['match']==1))
            west = func1(linecat['center'][ind])
            #west = dln.poly(linecat['center'][ind],coef1)            
            medoff = np.median(linecat['wave'][ind]-west)
            woff1[i] = medoff
            wavecorr1[ind] = linecat['wave'][ind]-medoff
        print('Offsets ',woff1)

        # 3) Second polynomial fit to all lines in all fibers with initial offsets applied
        func2 = fitter(func_init,linecat['center'][gline],wavecorr1[gline])
        dwave2 = func2(linecat['center'][gline]) - wavecorr1[gline]
        #coef2 = dln.poly_fit(linecat['center'][gline],wavecorr1[gline],4)
        #dwave2 = dln.poly(linecat['center'][gline],coef2) - wavecorr1[gline]
        sigwave2 = dln.mad(dwave2)
        #print('Second fit coefficients ',func2.parameters)
        #print('Second fit coefficients ',coef2)        
        print('Second sigma = %f6.2' % sigwave2)

        # 4) outlier rejection
        # take sections and robust fit a line, then reject outliers
        xx = linecat['center'][gline]
        yy = wavecorr1[gline]
        si = np.argsort(xx)
        step = 256
        bad = np.zeros(len(gline),bool)
        for i in range(8):
            g, = np.where((xx>=i*step) & (xx<(i+1)*step))
            cc,absdev = dln.ladfit(xx[g],yy[g])
            cc = cc[::-1]  # flip
            diff = yy[g]-dln.poly(xx[g],cc)
            sig = dln.mad(diff)
            bd, = np.where(np.abs(diff) > 3*sig)
            if len(bd)>0:
                bad[g[bd]] = True
            #print(i,i*step,(i+1)*step,len(bd))
        bd, = np.where(bad==True)
        if len(bd)>0:
            print('Rejecting '+str(len(bd))+' outlier points')
            linecat['outlier'][gline[bd]] = True
            # Select new good points
            gline, = np.where((linecat['match']==1) & (linecat['outlier']==False))
            xx = linecat['center'][gline]
            yy = wavecorr1[gline]
            si = np.argsort(xx)
        
        # 5) Fit polynomial again and linear fits for each fiber
        func2 = fitter(func_init,xx[si],yy[si])
        #coef2 = dln.poly_fit(xx[si],yy[si],4)      
        woff2 = np.zeros((2,nap),float)
        wavecorr2 = np.zeros(nlinecat,float)
        for i,ap in enumerate(gdap):
            ind, = np.where((linecat['aperture']==ap) & (linecat['match']==1))
            wdiff = linecat['wave'][ind] - func2(linecat['center'][ind])         
            #wdiff = linecat['wave'][ind] - dln.poly(linecat['center'][ind],coef2)
            wcc1,wabsdev1 = dln.ladfit(linecat['center'][ind],wdiff)
            wcc1 = wcc1[::-1]  # flip
            # Do outlier rejection
            wsig = dln.mad(wdiff-dln.poly(linecat['center'][ind],wcc1))
            wgd, = np.where(np.abs(wdiff-dln.poly(linecat['center'][ind],wcc1)) < 4*wsig)
            print(i,len(ind)-len(wgd))
            wcc,wabsdev = dln.ladfit(linecat['center'][ind][wgd],wdiff[wgd])
            wcc = wcc[::-1]  # flip
            wavecorr2[ind] = linecat['wave'][ind] - dln.poly(linecat['center'][ind],wcc)
            # Convert to Legendre coefficents by fitting to the poly fit
            func_init1 = models.Legendre1D(degree=1,domain=[0,npix])
            fitter1 = fitting.LinearLSQFitter()
            wfunc1 = fitter1(func_init1,x,dln.poly(x,wcc))
            woff2[:,i] = wfunc1.parameters
            #if diag==True:
            #    plt.scatter(linecat['center'][ind],wdiff)
            #    plt.scatter(linecat['center'][ind][wgd],wdiff[wgd])            
            #    plt.plot(linecat['center'][ind],dln.poly(linecat['center'][ind],wcc))
            #    plt.plot(linecat['center'][ind],wfunc1(linecat['center'][ind]))
            #    plt.show()
            
        # 6) Another outlier rejection
        yy = wavecorr2[gline]
        func3 = fitter(func_init,xx[si],yy[si])
        dwave2 = yy - func3(xx)
        #dwave2 = yy - dln.poly(xx,coef2)        
        sig2 = dln.mad(dwave2)
        bd, = np.where(np.abs(dwave2) > 3*sig2)
        if len(bd)>0:
            print('Rejecting '+str(len(bd))+' more outlier points')
            linecat['outlier'][gline[bd]] = True
            # Select new good points
            gline, = np.where((linecat['match']==1) & (linecat['outlier']==False))
            xx = linecat['center'][gline]
            yy = wavecorr2[gline]
            si = np.argsort(xx)

        # 7) Fit polynomial again
        func4 = fitter(func_init,xx[si],yy[si])
        dwave4 = yy - func4(xx)
        #coef4 = dln.poly_fit(xx[si],yy[si],5)
        #dwave4 = yy - dln.poly(xx,coef4)
        sig4 = dln.mad(dwave4)
        print('sigma = '+str(sig4))

        # 8) Fit linear fits for each fiber one last time
        woff3 = np.zeros((2,nap),float)
        wavecorr3 = np.zeros(nlinecat,float)
        for i,ap in enumerate(gdap):
            ind, = np.where((linecat['aperture']==ap) & (linecat['match']==1))
            wdiff = linecat['wave'][ind] - func4(linecat['center'][ind])
            #wdiff = dwave4[ind]
            wcc1,wabsdev1 = dln.ladfit(linecat['center'][ind],wdiff)
            wcc1 = wcc1[::-1]  # flip
            # Do outlier rejection
            wsig = dln.mad(wdiff-dln.poly(linecat['center'][ind],wcc1))
            wgd, = np.where(np.abs(wdiff-dln.poly(linecat['center'][ind],wcc1)) < 4*wsig)
            print(i,len(ind)-len(wgd))
            wcc,wabsdev = dln.ladfit(linecat['center'][ind][wgd],wdiff[wgd])
            wcc = wcc[::-1]  # flip
            wavecorr3[ind] = linecat['wave'][ind] - dln.poly(linecat['center'][ind],wcc)
            # Convert to Legendre coefficents by fitting to the poly fit
            func_init1 = models.Legendre1D(degree=1,domain=[0,npix])
            fitter1 = fitting.LinearLSQFitter()
            wfunc1 = fitter1(func_init1,x,dln.poly(x,wcc))
            woff3[:,i] = wfunc1.parameters
            if diag==True:
                plt.scatter(linecat['center'][ind],wdiff)
                plt.scatter(linecat['center'][ind][wgd],wdiff[wgd])            
                plt.plot(linecat['center'][ind],dln.poly(linecat['center'][ind],wcc))
                plt.plot(linecat['center'][ind],wfunc1(linecat['center'][ind]))
                plt.xlabel('X')
                plt.ylabel('Wavelength Residuals')
                plt.title('Aperture ')
                plt.show()
        

        # Should we try to ID more lines??

        #if 61 in gdap:
        #    import pdb; pdb.set_trace()
        # aperture 61 has some issues
        # but ThAr spectrum is not great
        # might need to fit it with the sky lines later

        
        # B-spline fit
        # smaller s gives more knots, len(spl.get_knots())
        # the knots are not equally spaced        
        #spl = UnivariateSpline(xx[si],yy[si],k=2,s=6)
        #nknots = 5
        #dknots = npix/nknots
        #knots = np.arange(nknots)*dknots + 0.5*dknots
        #spl = LSQUnivariateSpline(xx[si],yy[si],knots,k=2)
        #dwave2 = spl(xx) - wavecorr1[gline]
        #xs = np.linspace(-3, 3, 1000)
        #plt.plot(xs, spl(xs), 'g', lw=3)

        yy = wavecorr3[gline]
        if diag==True:
            ##plt.plot(spl(x)-dln.poly(x,coef1))  
            ##plt.scatter(xx,dwave2,c=linecat['aperture'][gline])
            ##plt.scatter(xx,yy-dln.poly(xx,coef1),c=linecat['aperture'][gline])
            #plt.scatter(xx,yy-dln.poly(xx,coef3),c=linecat['aperture'][gline])
            plt.scatter(xx,yy-func4(xx),c=linecat['aperture'][gline])        
            plt.show()

        
        # Now get coefficients for each aperture
        func_init = models.Legendre1D(degree=order,domain=[0,npix])
        fitter = fitting.LinearLSQFitter()        
        for i,ap in enumerate(gdap):
            apind, = np.where(apinfo['aperture']==ap)
            idline1 = idlines[apind[0]]
            
            indall, = np.where(linecat['aperture']==ap)
            linecat1 = linecat[indall]
            ind, = np.where((linecat1['match']==1) & (linecat1['outlier']==False))

            #func = fitter(func_init,linecat1['center'][ind],linecat1['wave'][ind])
            #print(func.parameters)
            #dwave = func(linecat1['center'][ind]) - linecat1['wave'][ind]
            func = func4.copy()
            func.c0 += woff3[0,i]
            func.c1 += woff3[1,i]            
            ww = func(linecat1['center'][ind])
            rms = np.sqrt(np.mean((linecat1['wave'][ind]-ww)**2))
            npoints = len(ind)
            y = np.ma.masked_array(deepcopy(linecat1['wave'][ind]),mask=np.full((npoints),False,dtype=bool))
            print(ap,rms,npoints)

            if diag==True:
                plt.scatter(linecat1['center'][ind],linecat1['wave'][ind]-ww)
                plt.plot([0,2048],[0,0])
                plt.xlabel('X')
                plt.ylabel('Wavelength Residuals')
                plt.title('Aperture '+str(ap)+' RMS=%6.4f' % rms)
                plt.show()
            
            id_lines_pix = linecat1['center']
            id_lines_wav = linecat1['wave']
            id_lines_used = np.zeros(len(indall),bool)
            id_lines_used[ind] = True
            #order = [id_lines_template.func._order]
            #func,rms,npoints,y = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
            #                                       order,rejection_iterations,rejection_sigma)        

            # Make the output
            wav = np.ma.masked_array(np.full((len(idline1.fit_lines.fit)),-999,dtype='float'),
                                     mask=np.full((len(idline1.fit_lines.fit)),True,dtype=bool))
            for k in range(0,len(id_lines_pix)):
                if id_lines_used[k]==True:
                    wav[k] = id_lines_wav[k]
                    wav.mask[k] = False
                else:
                    wav.mask[k] = True
                    
            resolution_order = 1
            resolution_rejection_iterations = 10
            resolution,resolution_rms,resolution_npoints = m2fs.get_resolution(deepcopy(idline1.fit_lines),deepcopy(wav),
                                                                               resolution_order,resolution_rejection_iterations)
                
            newidline = m2fs.id_lines(aperture=idline1.aperture,fit_lines=idline1.fit_lines,wav=wav,func=func,
                                      rms=rms,npoints=npoints,resolution=resolution,
                                      resolution_rms=resolution_rms,resolution_npoints=resolution_npoints)

            idlines[apind[0]] = newidline

        #import pdb; pdb.set_trace()
        
    if outfile=='':
        outfile = datadir+'/'+root0+'_blockwavefit.pickle'
    if os.path.exists(outfile): os.remove(outfile)
    print('Writing results to '+outfile)
    pickle.dump(idlines,open(outfile,'wb'))

    #import pdb; pdb.set_trace()
