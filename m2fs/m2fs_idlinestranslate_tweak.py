#!/usr/bin/env python

import os
import numpy as np
import dill as pickle
from astropy.io import fits
from dlnpyutils import bindata, utils as dln
from argparse import ArgumentParser
from doppler.rv import ccorrelate

# take flux, subtract scaled sky spectrum, then fit a continuum and subtract it
# then check the median (subtracted) flux of the pixels with bright sky lines
# want to minimize that

import astropy.units as u
from astropy.modeling import models,fitting
from copy import deepcopy
import copy
import m2fs_process as m2fs


if __name__ == "__main__":
    parser = ArgumentParser(description='Run Doppler fitting on spectra')
    parser.add_argument('idlinefile', type=str, nargs=1, help='M2FS ID lines translate filename')
    parser.add_argument('--outfile', type=str, nargs=1, default='', help='Output file')
    parser.add_argument('--directory', type=str, nargs=1, default='./', help='Main directory with the data')    
    args = parser.parse_args()
    idlinefile = args.idlinefile[0]
    outfile = dln.first_el(args.outfile)
    directory = dln.first_el(args.directory)

    npix = 2048
    x = np.arange(npix)
    
    print('Fixing wavelength solution for '+idlinefile)

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
    aperture = np.zeros(n,int)
    rms = np.zeros(n,float)
    npoints = np.zeros(n,int)
    for i in range(n):
        aperture[i] = idlines[i].aperture
        rms[i] = idlines[i].rms
        npoints[i] = idlines[i].npoints        
    used, = np.where(npoints>0)
    medrms = np.median(rms[used])
    mednpoints = np.median(npoints[used])

    # Bad ones to fix
    #badap1, = np.where((npoints>0) & ((npoints<0.5*mednpoints) | (rms>3*medrms)))
    goodap1, = np.where((npoints>0.5*mednpoints) & (rms<1.5*medrms))

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
    goodap, = np.where((npoints>0) & (sig<1.0) & (std<1.0) & (rms<0.1) & (npoints>20))
    badap, = np.where((npoints>0) & ((sig>=1.0) | (std>=1.0) | (rms>=0.1) | (npoints<=20)))    
    nbad = len(badap)
    print(str(nbad)+' bad apertures to fix')

    #goodap, = np.where((npoints>0) & (sig<1.0) & (std<1.0) & (npoints>0.5*mednpoints) & (rms<0.08))
    #badap, = np.where(rms > 0.08)
    #nbad = len(badap)
    #import pdb; pdb.set_trace()
    
    # Loop over bad apertures
    origidlines = idlines.copy()
    for i in range(nbad):
        idline1 = idlines[badap[i]]
        print(' ')
        print('Fixing aperture '+str(idline1.aperture))
        spec = extract1d_array[badap[i]]
        # Find a neighboring good aperture
        bestind = goodap[np.argmin(np.abs(aperture[badap[i]]-aperture[goodap]))]
        #if idline1.aperture==61:
        #    bestind = 62
        print('Best reference aperture: '+str(aperture[bestind]))
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

        # Perform an initial fit
        gline1, = np.where(linecat['match']==1)
        meanwave = np.mean(linecat['wave'][gline1])
        coef1 = dln.poly_fit(linecat['center'][gline1],linecat['wave'][gline1],3)
        dwave = dln.poly(linecat['center'][gline1],coef1) - linecat['wave'][gline1]
        sigwave = dln.mad(dwave)
        if sigwave>0.1:
            print('Sigma high.  Trying higher order')
            coef1 = dln.poly_fit(linecat['center'][gline1],linecat['wave'][gline1],4)
            dwave = dln.poly(linecat['center'][gline1],coef1) - linecat['wave'][gline1]
            sigwave = dln.mad(dwave)
        print('Sigma of initial fit %6.4f A' % sigwave)
        bdline, = np.where(np.abs(dwave) > 3*sigwave)
        if len(bdline)>0:
            linecat['outlier'][gline1[bdline]] = True
            print('throwing out '+str(len(bdline))+' outlier lines')

        #if sigwave>0.1:
        #    print('Stopping. sigma too high')
        #    import pdb; pdb.set_trace()
            
        # Now fit the function
        gline, = np.where((linecat['match']==1) & (linecat['outlier']==False))
        id_lines_pix = linecat['center']
        id_lines_wav = linecat['wave']
        id_lines_used = np.zeros(nlines,bool)
        id_lines_used[gline] = True

        # Try many different orders
        rejection_iterations = [10]
        rejection_sigma = [3.0]
        #order = [id_lines_template.func._order]
        order = [3]
        func3,rms3,npoints3,y3 = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
                                                   order,rejection_iterations,rejection_sigma)
        order = [4]
        func4,rms4,npoints4,y4 = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
                                                   order,rejection_iterations,rejection_sigma)
        order = [5]
        func5,rms5,npoints5,y5 = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
                                                   order,rejection_iterations,rejection_sigma)        
        # Pick best one
        funcarr = [func3,func4,func5]
        rmsarr = [rms3,rms4,rms5]
        npointsarr = [npoints3,npoints4,npoints5]
        orderarr = [3,4,5]
        yarr = [y3,y4,y5]
        bestind = np.argmin(np.array(rmsarr))
        order = [ orderarr[bestind] ]
        func = funcarr[bestind]
        rms = rmsarr[bestind]
        npoints = npointsarr[bestind]
        y = yarr[bestind]
        print('best order = '+str(order[0]))
        print('rms = %6.3f' % rms)

        if rms>0.1:
            print('Stopping. rms too high')
            import pdb; pdb.set_trace()
        
        
        ## Lower order if necessary
        #if rms>0.1:
        #    print('RMS too high.  Lowering order')
        #    order = [order[0]-1]
        #    func,rms,npoints,y = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
        #                                           order,rejection_iterations,rejection_sigma)

        ## Lower order if necessary
        #if rms>0.1:
        #    print('RMS too high.  Lowering order')
        #    order = [order[0]-1]
        #    func,rms,npoints,y = m2fs.id_lines_fit(id_lines_pix[gline],id_lines_wav[gline],id_lines_used[gline],
        #                                           order,rejection_iterations,rejection_sigma)

            
        # Make the output
        wav = np.ma.masked_array(np.full((len(idline1.fit_lines.fit)),-999,dtype='float'),mask=np.full((len(idline1.fit_lines.fit)),True,dtype=bool))
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
        idlines[badap[i]] = newidline

        #import pdb; pdb.set_trace()
        

    #outfile = 'test.pickle'
    if outfile=='':
        outfile = datadir+'/'+root0+'_id_lines_array1.pickle'
    if os.path.exists(outfile): os.remove(outfile)
    print('Writing results to '+outfile)
    pickle.dump(idlines,open(outfile,'wb'))

    import pdb; pdb.set_trace()
