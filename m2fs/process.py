#!/usr/bin/env python

def get_datadir(m2fsrun):
    maindir = '/nfs/nas-0-9/mgwalker.proj/m2fs/m2fs.astro.lsa.umich.edu/data/'
    if m2fsrun=='feb14': datadir = maindir+'Feb2014/'
    if m2fsrun=='dec14': datadir = maindir+'Dec2014/'
    if m2fsrun=='feb15': datadir = maindir+'Feb2015/'
    if m2fsrun=='jul15': datadir = maindir+'Jul2015/'
    if m2fsrun=='sep15': datadir = maindir+'Sep2015/'
    if m2fsrun=='nov15': datadir = maindir+'Nov2015/'
    if m2fsrun=='feb16': datadir = maindir+'Feb2016/'
    if m2fsrun=='jun16': datadir = maindir+'Jun2016/'
    if m2fsrun=='aug16': datadir = maindir+'AugSep2016/'
    if m2fsrun=='nov16': datadir = maindir+'NovDec2016/'
    if m2fsrun=='feb17': datadir = maindir+'FebMar2017/'
    if m2fsrun=='may17': datadir = maindir+'MayJun2017/'
    if m2fsrun=='sep17': datadir = maindir+'Sep2017/'
    if m2fsrun=='nov17': datadir = maindir+'Nov2017/'
    if m2fsrun=='feb18': datadir = maindir+'Feb2018/'
    if m2fsrun=='may18': datadir = maindir+'May2018/'
    if m2fsrun=='aug18': datadir = maindir+'Aug2018/'
    if m2fsrun=='nov18': datadir = maindir+'NovDec2018/'
    if m2fsrun=='feb19': datadir = maindir+'FebMar2019/'
    if m2fsrun=='may19': datadir = maindir+'MayJun2019/'
    if m2fsrun=='aug19': datadir = maindir+'AugSep2019/'
    if m2fsrun=='nov19': datadir = maindir+'Nov2019/'
    if m2fsrun=='jan20': datadir = maindir+'JanFeb2020/'
    if m2fsrun=='bug16': datadir = maindir+'AugSep2016/'
    return(datadir)

def get_resolution(fit_lines,wav,resolution_order,resolution_rejection_iterations):
    import numpy as np
    import scipy
    import astropy.units as u
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    import matplotlib.pyplot as plt

    resolution_init = models.Polynomial1D(degree=resolution_order)
    resolution_fitter = fitting.LinearLSQFitter()
    x = np.array([fit_lines.fit[q].mean.value for q in range(0,len(fit_lines.fit))])
    resolution_y = np.ma.masked_array([fit_lines.fit[q].stddev.value for q in range(0,len(fit_lines.fit))],mask=wav.mask)
    resolution = resolution_init
    use = np.where(resolution_y.mask==False)[0]
    resolution_rms = np.sqrt(np.mean((resolution_y[use]-resolution(x[use]))**2))
    for q in range(0,resolution_rejection_iterations):
        use = np.where(resolution_y.mask==False)[0]
        if len(use)>3.:
            resolution = resolution_fitter(resolution_init,x[use],resolution_y[use])
            resolution_rms = np.sqrt(np.mean((resolution_y[use]-resolution(x[use]))**2))
            outlier = (np.where(np.abs(resolution_y-resolution(x))>3.*resolution_rms))[0]  # reject outliers
            resolution_y.mask[outlier] = True
    use = np.where(resolution_y.mask==False)[0]
    resolution_npoints = len(use)

#    plt.scatter(x,resolution_y,s=3,color='k')
#    plt.scatter(x[use],resolution_y[use],s=3,color='g')
#    plt.plot(x,resolution(x),color='r')
#    plt.xlim([np.min(x[use]),np.max(x[use])])
#    plt.show()
#    plt.close()
    return resolution,resolution_rms,resolution_npoints

def line_id_add_lines(linelist,line_centers,id_lines_used,func0,id_lines_tol_angs):
    import numpy as np
    new_pix = []
    new_wav = []
    new_used = []
    for j in range(0,len(line_centers)):
        if j not in id_lines_used:
            dist = np.sqrt((linelist.wavelength-func0(line_centers[j]))**2)
            nan = np.where(dist!=dist)[0]
            if len(nan)>0:
                dist[nan] = np.inf  # artificially replace any nans with huge values
            best = np.where(dist==np.min(dist))[0][0]
            if dist[best]<id_lines_tol_angs:
                print(line_centers[j])
                new_pix.append(line_centers[j])
                new_wav.append(np.float(linelist.wavelength[best]))
                new_used.append(j)
    print('added ',len(new_used),' new lines from line list')

    return new_pix,new_wav,new_used

def plot_trace(j,trace_x,trace_y,trace_z,y,trace_func,pix_min,pix_max,rms,npoints,order,fig):
    ax1 = fig.add_subplot(211)
    ax1.cla()
    ax1.scatter(trace_x,trace_y,s=3,color='k')
    ax1.scatter(trace_x[y.mask==True],trace_y[y.mask==True],s=13,marker='x',color='y')
    ax1.axvline(x=pix_min,color='k',lw=0.3,linestyle=':')
    ax1.axvline(x=pix_max,color='k',lw=0.3,linestyle=':')
    ax1.set_xlabel('trace x [pixel]')
    ax1.set_ylabel('trace y [pixel]')
    ax1.plot(trace_x,trace_func(trace_x),color='r')
    ax1.text(0,1.01,'n_points='+str(npoints),transform=ax1.transAxes)
    ax1.text(0,1.06,'order='+str(order),transform=ax1.transAxes)
    ax1.text(0,1.11,'rms='+str(rms),transform=ax1.transAxes)
    ax1.text(1,1.01,'aperture='+str(j+1),horizontalalignment='right',transform=ax1.transAxes)

    ax2 = fig.add_subplot(212)
    ax2.cla()
    ax2.scatter(trace_x,trace_z,s=3,color='k')
    ax2.scatter(trace_x[y.mask==True],trace_z[y.mask==True],s=13,marker='x',color='y')
    ax2.axvline(x=pix_min,color='k',lw=0.3,linestyle=':')
    ax2.axvline(x=pix_max,color='k',lw=0.3,linestyle=':')
    ax2.set_xlabel('trace x [pixel]')
    ax2.set_ylabel('counts')
    return ax1,ax2

def plot_id_lines(extract1d,continuum,fit_lines,line_centers,id_lines_pix,id_lines_wav,id_lines_used,func,fig):
    import numpy as np
    from astropy import units as u

    xlim = [np.min(extract1d.spec1d_pixel[extract1d.spec1d_mask==False]),np.max(extract1d.spec1d_pixel[extract1d.spec1d_mask==False])]
    ax1 = fig.add_subplot(211)
    ax1.cla()
    ax1.set_xlim([np.min(extract1d.spec1d_pixel[extract1d.spec1d_mask==False]),
                  np.max(extract1d.spec1d_pixel[extract1d.spec1d_mask==False])])
    ax1.plot(extract1d.spec1d_pixel[extract1d.spec1d_mask==False],
             extract1d.spec1d_flux[extract1d.spec1d_mask==False],color='k',lw=0.5)
    ax1.plot(extract1d.spec1d_pixel[extract1d.spec1d_mask==False],
             continuum(extract1d.spec1d_pixel[extract1d.spec1d_mask==False]),color='b',lw=0.3)
    ax1.set_xlabel('pixel')
    ax1.set_ylabel('counts')
    ax1.set_xscale('linear')
    ax1.set_yscale('linear')
    ax1.set_xlim(xlim)
    for j in range(0,len(fit_lines.fit)):
        use = np.where((extract1d.spec1d_pixel>=fit_lines.subregion[j].lower.value)&(extract1d.spec1d_pixel<=fit_lines.subregion[j].upper.value))[0]
        sub_spectrum_pixel = extract1d.spec1d_pixel[use]*u.AA
        y_fit = fit_lines.fit[j](sub_spectrum_pixel).value+continuum(sub_spectrum_pixel.value)
        ax1.plot(sub_spectrum_pixel,y_fit,color='r',lw=0.3)
        ax1.axvline(x=fit_lines.fit[j].mean.value,linestyle=':',color='k',lw=0.3)
    for j in range(0,len(id_lines_pix)):
         ax1.axvline(x=id_lines_pix[j],color='g',lw=1,linestyle='-')
#    ax1.text(0,1.01,'n_points='+str(npoints),transform=ax1.transAxes)
#    ax1.text(0,1.06,'order='+str(order),transform=ax1.transAxes)
#    ax1.text(0,1.11,'rms='+str(rms),transform=ax1.transAxes)
#    ax1.text(1,1.01,'aperture='+str(j+1),horizontalalignment='right',transform=ax1.transAxes)
    ax1.set_xscale('linear')
    ax1.set_yscale('linear')

    ax2 = fig.add_subplot(212)
    ax2.cla()
#    ax2.scatter(id_lines_pix,id_lines_wav,color='k',s=3)
    ax2.scatter(id_lines_pix,id_lines_wav-func(id_lines_pix),color='k',s=3)
    x = extract1d.spec1d_pixel[extract1d.spec1d_mask==False]
#    ax2.plot(x,func(x),color='r')
    ax2.set_xlim(xlim)
    ax2.set_xlabel('pixel')
    ax2.set_ylabel(r'$\lambda$ [Angstroms]')

    return ax1,ax2

def id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,order,rejection_iterations,rejection_sigma):
    import numpy as np
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    from copy import deepcopy
#    print('fitting wavelength solution')

    x = np.array(id_lines_pix)
    y = np.ma.masked_array(deepcopy(id_lines_wav),mask=np.full((len(id_lines_wav)),False,dtype=bool))
    for q in y.mask:
        y.mask[q] = False
    order0 = np.max([1,np.min([order[len(order)-1],len(x)-2])])
    func_init = models.Legendre1D(degree=order0)
    fitter = fitting.LinearLSQFitter()
    for q in range(0,rejection_iterations[len(rejection_iterations)-1]):
        use = np.where(y.mask==False)[0]
        order0 = np.max([1,np.min([order[len(order)-1],len(use)-2])])
        func_init = models.Legendre1D(degree=order0)
#        print(id_lines_wav[use],' use')
        if len(use)>1:
            func0 = fitter(func_init,x[use],y[use])
        else:
            func0 = func_init
        rms0 = np.sqrt(np.mean((y[use]-func0(x)[use])**2))
        outlier = np.where(np.abs(y-func0(x))>rejection_sigma[len(rejection_sigma)-1]*rms0)[0]
        y.mask[outlier] = True
    use = np.where(y.mask==False)[0]
    npoints0 = len(use)
#        print('order of polynomial fit is ',np.min([order[len(order)-1],len(use)]))
    print('order of polynomial fit is ',order0)
    print('number of points in fit=',npoints0)
    print('rejected points = ',len(np.where(y.mask==True)[0]))
    print('rms='+str(rms0))
#    print(func0)
    print('')
    return  func0,rms0,npoints0,y

def plot_boundary(lower_x,lower_y,upper_x,upper_y,reject_x,reject_y,data,fig):
    import numpy as np
    def getxy(x0,y0):
        x = x0
        y = y0
        x = np.array(x)
        y = np.array(y)
        order = np.argsort(y)
        x = x[order]
        y = y[order]
        x = np.append(x[0],x)
        y = np.append(0,y)
        x = np.append(x,x[len(x)-1])
        y = np.append(y,len(data.data))
        return x,y
    ax1 = fig.add_subplot(111)
    ax1.cla()

    vmin = np.quantile(data.data.flatten(),0.05)
    vmax = np.quantile(data.data.flatten(),0.95)
    ax1.imshow(data.data,cmap='gray',vmin=vmin,vmax=vmax)
    if len(lower_x)>0:
        print(lower_x)
        print(lower_y)
        x,y = getxy(lower_x,lower_y)
        print(x,y)
        ax1.scatter(x,y,color='y',s=10,alpha=0.2)
        ax1.plot(x,y,color='y',alpha=0.2)
    if len(upper_x)>0:
        x,y = getxy(upper_x,upper_y)
        ax1.scatter(x,y,color='y',s=10,alpha=0.2)
        ax1.plot(x,y,color='y',alpha=0.2)
    if len(reject_x)>0:
        x,y = getxy(reject_x,reject_y)
        ax1.scatter(x,y,color='r',s=10,alpha=0.2)
        ax1.plot(x,y,color='r',alpha=0.2)

def on_key_boundary(event,args_list):
    import numpy as np
    #import m2fs_process as m2fs
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D

    global lower_x,lower_y,upper_x,upper_y,reject_x,reject_y
    data,lower_x,lower_y,upper_x,upper_y,reject_x,reject_y,fig=args_list
    if event.key=='l':  # add lower boundary point
        print('add boundary point at (',event.xdata,event.ydata,')')
        lower_x.append(event.xdata)
        lower_y.append(event.ydata)

    if event.key=='u':  # add upper boundary point
        print('add boundary point at (',event.xdata,event.ydata,')')
        upper_x.append(event.xdata)
        upper_y.append(event.ydata)

    if event.key=='b':  # add upper boundary point
        print('reject point at (',event.xdata,event.ydata,')')
        reject_x.append(event.xdata)
        reject_y.append(event.ydata)

    if event.key=='d':  # delete nearest boundary point
        print('delete point nearest (',event.xdata,event.ydata,')')
        x0 = event.xdata
        y0 = event.ydata
        if len(lower_x)>0:
            dist_lower = (lower_x-x0)**2+(lower_y-y0)**2
            min1 = np.min(dist_lower)
            best1 = np.where(dist_lower==min1)[0][0]
        else:
            min1 = 1.e+30
        if len(upper_x)>0:
            dist_upper = (upper_x-x0)**2+(upper_y-y0)**2
            min2 = np.min(dist_upper)
            best2 = np.where(dist_upper==min2)[0][0]
        else:
            min2 = 1.e+30
        if min1<=min2:
            if len(lower_x)>0:
                del lower_x[best1]
                del lower_y[best1]
        else:
            if len(upper_x)>0:
                del upper_x[best2]
                del upper_y[best2]

    ax1 = plot_boundary(lower_x,lower_y,upper_x,upper_y,reject_x,reject_y,data,fig)
    fig.canvas.draw_idle()
    return

def on_key_find(event,args_list):
    import numpy as np
    #import m2fs_process as m2fs
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D

    print('you pressed ',event.key)

    global columnspec_array,subregion,fit,realvirtual,initial,window
    columnspec_array,column,subregion,fit,realvirtual,initial,window,fig = args_list

    if event.key=='z':  # revert to initial apertures
        to_delete = np.where(np.logical_not(initial))[0]
        j = 0
        for i in range(0,len(to_delete)):
            del initial[to_delete[i]-j]
            del subregion[to_delete[i]-j]
            del fit[to_delete[i]-j]
            del realvirtual[to_delete[i]-j]
            j += 1
        subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)

    if event.key=='d':
        aperture_centers = [fit[q].mean for q in range(0,len(fit))]
        x0 = event.xdata
        y0 = event.ydata
        dist = (aperture_centers-event.xdata)**2
        best = np.where(dist==np.min(dist))[0][0]
        print('deleting aperture '+str(best+1)+', centered at '+str(aperture_centers[best]))
        del subregion[best]
        del fit[best]
        del realvirtual[best]
        del initial[best]
        subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)

    if event.key=='e':
        aperture_centers = [fit[q].mean for q in range(0,len(fit))]
        print('deleting all apertures ')
        for q in range(1,len(subregion)):
            del(subregion[len(subregion)-1])
        for q in range(1,len(fit)):
            del(fit[len(fit)-1])
        for q in range(1,len(realvirtual)):
            del(realvirtual[len(realvirtual)-1])
        for q in range(1,len(initial)):
            del(initial[len(initial)-1])
        subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)

    if event.key=='n':
        new_center = np.float(event.xdata)
        x_center = new_center
        spec1d = Spectrum1D(spectral_axis=columnspec_array[column].pixel,
                            flux=columnspec_array[column].spec*u.electron,
                            uncertainty=columnspec_array[column].err,
                            mask=columnspec_array[column].mask)
        subregion0,fit0 = fit_aperture(spec1d-columnspec_array[column].continuum(columnspec_array[column].pixel.value),
                                       window,x_center)
        subregion.append(subregion0)
        fit.append(fit0)
        realvirtual.append(True)
        initial.append(False)
        subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)

    if event.key=='a':
        new_center = np.float(event.xdata)
        x_center = new_center
        val1 = x_center-window/2.
        val2 = x_center+window/2.
        subregion.append(SpectralRegion(val1*u.AA,val2*u.AA))  # define extraction region from window
        aaa = np.float(np.max(columnspec_array[column].spec-columnspec_array[column].continuum(columnspec_array[column].pixel.value)))
        halfwindow = window/2.
        fit.append(models.Gaussian1D(amplitude=aaa*u.electron,
                                     mean=x_center*u.AA,stddev=halfwindow*u.AA))
        realvirtual.append(False)
        initial.append(False)
        subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)

    x = np.arange(0,len(columnspec_array[column].spec))
    ax1 = fig.add_subplot(111)
    ax1.cla()
    ax1.plot(columnspec_array[column].pixel,columnspec_array[column].spec,lw=0.3,color='k')
    ax1.plot(columnspec_array[column].pixel,columnspec_array[column].continuum(x),color='green',lw=1)
    ax1.set_xlim([0,np.max(x)])
    for j in range(0,len(fit)):
        use = np.where((columnspec_array[column].pixel.value>=subregion[j].lower.value)&(columnspec_array[column].pixel.value<=subregion[j].upper.value))[0]
        sub_spectrum_pixel = columnspec_array[column].pixel[use]
        y_fit = fit[j](sub_spectrum_pixel).value+columnspec_array[column].continuum(sub_spectrum_pixel.value)

        ax1.axvline(x=fit[j].mean.value,color='b',lw=0.3,alpha=0.7)
        ax1.plot(sub_spectrum_pixel.value,y_fit,color='r',lw=0.2,alpha=0.7)
        ax1.text(fit[j].mean.value,0,str(j+1),fontsize=8)
        ax1.set_ylim([0,10000])
        if not realvirtual[j]:
            ax1.axvline(x=fit[j].mean.value,color='k',lw=0.5,alpha=1,linestyle='--')
    fig.canvas.draw_idle()
    return

def on_key_trace(event,args_list):
    import numpy as np
    #import m2fs_process as m2fs
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    print('you pressed ',event.key)

    global trace_x,trace_y,trace_z,y,sig_y,pix_min,pix_max,func,rms,npoints,order,rejection_sigma,trace_rejection_iterations,trace_func_init
    j,trace_x,trace_y,trace_z,y,sig_y,pix_min,pix_max,func,rms,npoints,order,rejection_sigma,trace_rejection_iterations,trace_func_init,fig = args_list

    if event.key=='d':
        x0 = event.xdata
        y0 = event.ydata
        candidates = np.where(y.mask==False)[0]
        dist = ((trace_x[candidates]-event.xdata)**2+(trace_y[candidates]-event.ydata)**2)
        best = np.where(dist==np.min(dist))[0][0]
        print('deleting point centered at ('+str(trace_x[candidates[best]])+','+str(trace_y[candidates[best]])+')')
        y.mask[candidates[best]] = True

    if event.key=='c':
        x0 = event.xdata
        y0 = event.ydata
        candidates = np.where(y.mask==True)[0]
        dist = ((trace_x[candidates]-event.xdata)**2+(trace_y[candidates]-event.ydata)**2)
        best = np.where(dist==np.min(dist))[0][0]
        print('re-instating point '+str(best+1)+', centered at ('+str(trace_x[best])+','+str(trace_y[best])+')')
        y.mask[candidates[best]] = False

    if event.key=='l':
        pix_min.append(event.xdata)
    if event.key=='u':
        pix_max.append(event.xdata)
    if event.key=='o':
        print('order of polynomial fit is ',order[len(order)-1])
        command = input('enter new order (must be integer): ')
        order.append(np.long(command))
    if event.key=='r':
        print('rejection sigma is ',rejection_sigma[len(rejection_sigma)-1])
        command = input('enter new rejection sigma (float): ')
        rejection_sigma.append(np.float(command))
    if event.key=='z':
        for q in range(1,len(pix_min)):
            del(pix_min[len(pix_min)-1])
        for q in range(1,len(pix_max)):
            del(pix_max[len(pix_max)-1])
        for q in range(1,len(order)):
            del(order[len(order)-1])
        for q in range(1,len(rejection_sigma)):
            del(rejection_sigma[len(rejection_sigma)-1])
        for q in range(0,len(y)):
            y.mask[q] = False

    y.mask[trace_x<pix_min[len(pix_min)-1]] = True
    y.mask[trace_x>pix_max[len(pix_max)-1]] = True
    fitter = fitting.LinearLSQFitter()
    for q in range(0,trace_rejection_iterations):
        use = np.where(y.mask==False)[0]
        trace_func = fitter(trace_func_init,trace_x[use],trace_y[use])
        trace_rms = np.sqrt(np.mean((y[use]-trace_func(trace_x)[use])**2))
        outlier = np.where(np.abs(trace_y-trace_func(trace_x))>rejection_sigma[len(rejection_sigma)-1]*trace_rms)[0]
        y.mask[outlier] = True

    use = np.where(y.mask==False)[0]
    trace_npoints = len(use)

    func.append(trace_func)
    rms.append(trace_rms)
    npoints.append(trace_npoints)
    print('fitting range: ',pix_min[len(pix_min)-1],pix_max[len(pix_max)-1])
    print('order of polynomial fit is ',order[len(order)-1])
    print('number of points in fit=',trace_npoints)
    print('rms='+str(trace_rms))
    print('')
    ax1,ax2 = plot_trace(j,trace_x,trace_y,trace_z,y,trace_func,pix_min[len(pix_min)-1],pix_max[len(pix_max)-1],
                         trace_rms,trace_npoints,order[len(order)-1],fig)
    fig.canvas.draw_idle()
    return

def on_key_id_lines(event,args_list):
    import numpy as np
    #import m2fs_process as m2fs
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    from copy import deepcopy

    print('you pressed ',event.key)

    global id_lines_pix,id_lines_wav,id_lines_used,order,rejection_iterations,rejection_sigma,func,rms,npoints
    extract1d,continuum,fit_lines,linelist,line_centers,id_lines_pix,id_lines_wav,id_lines_used,order,rejection_iterations,rejection_sigma,func,rms,npoints,id_lines_tol_angs,fig = args_list

    if event.key=='.':
        x0 = event.xdata
        y0 = event.ydata
        print('cursor pixel = ',str(x0))
        dist = (line_centers-x0)**2
        best = np.where(dist==np.min(dist))[0][0]
        print('nearest line pixel = ',str(line_centers[best]))
        used = np.array(np.where(id_lines_pix==line_centers[best])[0])
        if len(used)>0:
            print('wavelength from line list =',str(np.array(id_lines_wav)[used][0]))
        else:
            print('line not yet identified')
        print('function wavelength =',str(func[len(func)-1](line_centers[best])))

    if event.key=='m':
        x0 = event.xdata
        y0 = event.ydata
        print('cursor pixel = ',str(x0))
        dist = (line_centers-x0)**2
        best = np.where(dist==np.min(dist))[0][0]
        print('nearest line pixel = ',str(line_centers[best]))
        command = input('enter wavelength in Angstroms \n')
        if command=='':
            print('no information entered')
        else:
            id_lines_pix.append(line_centers[best])
            id_lines_wav.append(np.float(command))
            id_lines_used.append(best)

    if event.key=='d':  # delete nearest boundary point
        print('delete point nearest (',event.xdata,event.ydata,')')
        if len(id_lines_pix)>0:
            dist = (id_lines_pix-event.xdata)**2
            best = np.where(dist==np.min(dist))[0][0]
            del id_lines_pix[best]
            del id_lines_wav[best]
            del id_lines_used[best]
        else:
            print('no ID\'d lines to delete!')

    if event.key=='o':
        print('order of polynomial fit is ',order[len(order)-1])
        command = input('enter new order (must be integer): ')
        if command=='':
            print('keeping original value')
        else:
            order.append(np.long(command))
        func0,rms0,npoints0,y = id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,order,
                                             rejection_iterations,rejection_sigma)
        func.append(func0)
        rms.append(rms0)
        npoints.append(npoints0)

    if event.key=='r':
        print('rejection sigma is ',rejection_sigma[len(rejection_sigma)-1])
        command = input('enter new rejection sigma (float): ')
        if command=='':
            print('keeping original value')
        else:
            rejection_sigma.append(np.float(command))
        func0,rms0,npoints0,y = id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,order,
                                             rejection_iterations,rejection_sigma)
        func.append(func0)
        rms.append(rms0)
        npoints.append(npoints0)

    if event.key=='t':
        print('number of rejection iterations is ',rejection_iterations[len(rejection_iterations)-1])
        command = input('enter new rejection iterations = (integer): ')
        if command=='':
            print('keeping original value')
        else:
            rejection_iterations.append(np.long(command))
        func0,rms0,npoints0,y = id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,order,
                                             rejection_iterations,rejection_sigma)
        func.append(func0)
        rms.append(rms0)
        npoints.append(npoints0)

    if event.key=='z':
        for q in range(1,len(id_lines_pix)):
            del(id_lines_pix[len(id_lines_pix)-1])
        for q in range(1,len(id_lines_wav)):
            del(id_lines_wav[len(id_lines_wav)-1])
        for q in range(1,len(id_lines_used)):
            del(id_lines_used[len(id_lines_used)-1])
        for q in range(1,len(order)):
            del(order[len(order)-1])
        for q in range(1,len(rejection_iterations)):
            del(rejection_iterations[len(rejection_iterations)-1])
        for q in range(1,len(rejection_sigma)):
            del(rejection_sigma[len(rejection_sigma)-1])
        for q in range(1,len(func)):
            del(func[len(func)-1])

    if event.key=='g':
        func0,rms0,npoints0,y = id_lines_fit(deepcopy(id_lines_pix),deepcopy(id_lines_wav),
                                             deepcopy(id_lines_used),order,rejection_iterations,
                                             rejection_sigma)
        func.append(func0)
        rms.append(rms0)
        npoints.append(npoints0)

    if event.key=='l':
        func0,rms0,npoints0,y = id_lines_fit(deepcopy(id_lines_pix),deepcopy(id_lines_wav),
                                             deepcopy(id_lines_used),order,rejection_iterations,
                                             rejection_sigma)
        new_pix,new_wav,new_used = line_id_add_lines(deepcopy(linelist),deepcopy(line_centers),
                                                     deepcopy(id_lines_used),deepcopy(func0),
                                                     id_lines_tol_angs)

        for i in range(0,len(new_pix)):
            id_lines_pix.append(new_pix[i])
            id_lines_wav.append(new_wav[i])
            id_lines_used.append(new_used[i])
        func0,rms0,npoints0,y = id_lines_fit(deepcopy(id_lines_pix),deepcopy(id_lines_wav),
                                             deepcopy(id_lines_used),order,rejection_iterations,
                                             rejection_sigma)
        func.append(func0)
        rms.append(rms0)
        npoints.append(npoints0)

    ax1,ax2 = plot_id_lines(extract1d,continuum,fit_lines,line_centers,id_lines_pix,
                            id_lines_wav,id_lines_used,func[len(func)-1],fig)
    fig.canvas.draw_idle()
    return

def on_click(event):  # (event.xdata,event.ydata) give the user coordinates, event.x and event.y are internal coordinates
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single',event.button,
           event.x,event.y,event.xdata,event.ydata))
    click.append((event.xdata,event.ydata))

def get_continuum(spec1d,continuum_rejection_low,continuum_rejection_high,
                  continuum_rejection_iterations,continuum_rejection_order):
    import numpy as np
    import astropy.units as u
    from astropy.modeling import models,fitting
    from copy import deepcopy
    import copy
    import matplotlib.pyplot as plt

#    continuum_init=models.Polynomial1D(degree=10)
    continuum_init = models.Chebyshev1D(degree=continuum_rejection_order)
    fitter = fitting.LinearLSQFitter()
#    lamb=(np.arange(len(spec1d.data),dtype='float'))*u.AA#unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    lamb = spec1d.spectral_axis  # (np.arange(len(spec1d.spectral_axis),dtype='float'))*u.AA#unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    y = np.ma.masked_array(spec1d.flux,mask=spec1d.mask)
#    print(y.mask[y!=y])
#    print(np.max(y))
#    plt.plot(y)
#    plt.show()
#    plt.close()
    for q in range(0,continuum_rejection_iterations):
        if len(np.where(y.mask==False)[0])>0:
            continuum = fitter(continuum_init,lamb.value[y.mask==False],y[y.mask==False])
        else:
            continuum = continuum_init
#        continuum = fitter(continuum_init,lamb.value,y)
        rms = np.sqrt(np.mean((y.data.value[y.mask==False]-continuum(lamb.value)[y.mask==False])**2))    
        outlier = (np.where(spec1d.data-continuum(lamb.value)>continuum_rejection_high*rms))[0]  # reject only high outliers (emission lines)
        y.mask[outlier] = True
        outlier = (np.where(spec1d.data-continuum(lamb.value)<continuum_rejection_low*rms))[0]  # reject only high outliers (emission lines)
        y.mask[outlier] = True
#         outlier = (np.where(spec1d.data-continuum(lamb.value)>1.*spec1d.uncertainty.quantity.value))[0]#reject only high outliers (emission lines)
#        y.mask[outlier] = True
    return continuum,rms


def get_continuum_throughputcorr(spec1d,continuum_rejection_iterations):
    import numpy as np
    import astropy.units as u
    from astropy.modeling import models,fitting
    from copy import deepcopy
    import copy
    import matplotlib.pyplot as plt

#    continuum_init = models.Polynomial1D(degree=10)
    continuum_init = models.Chebyshev1D(degree=10)
    fitter = fitting.LinearLSQFitter()
#    lamb = (np.arange(len(spec1d.data),dtype='float'))*u.AA#unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    lamb = spec1d.spectral_axis  # (np.arange(len(spec1d.spectral_axis),dtype='float'))*u.AA#unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    y = np.ma.masked_array(spec1d.flux,mask=spec1d.mask)
    for q in range(0,continuum_rejection_iterations):
#        y = np.ma.masked_array(deepcopy(spec1d.data),mask=deepcopy(spec1d.mask))
        continuum = fitter(continuum_init,lamb.value[y.mask==False],y[y.mask==False])
        outlier = (np.where(np.abs(spec1d.data-continuum(lamb.value))>3.*spec1d.uncertainty.quantity.value))[0]  # reject only high outliers (emission lines)
        y.mask[outlier] = True
#        plt.scatter(lamb[y.mask==False],y[y.mask==False])
#        plt.plot(lamb.value,continuum(lamb.value))
#        plt.show()
#        plt.close()
    rms = np.sqrt(np.mean((y.data.value[y.mask==False]-continuum(lamb.value)[y.mask==False])**2))    
    return continuum,rms


def get_cr_reject(spec1d,cr_rejection_low,cr_rejection_high,cr_rejection_order,cr_rejection_iterations,
                  cr_rejection_tol,cr_rejection_collateral):
    import numpy as np
    import astropy.units as u
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    from copy import deepcopy
    import copy
    import matplotlib.pyplot as plt

    spec1d0 = Spectrum1D(spectral_axis=deepcopy(spec1d.spec1d_pixel)*u.AA,flux=deepcopy(spec1d.spec1d_flux),
                         uncertainty=deepcopy(spec1d.spec1d_uncertainty),mask=deepcopy(spec1d.spec1d_mask))
    y = np.ma.masked_array(spec1d.spec1d_flux.value,mask=spec1d.spec1d_mask)
    if len(np.where(spec1d0.mask==False)[0])>100.:
        continuum0,rms0 = get_continuum(spec1d0,cr_rejection_low,cr_rejection_high,
                                        cr_rejection_iterations,cr_rejection_order)
 
#    plt.plot(spec1d.spec1d_pixel[spec1d.spec1d_mask==False],spec1d.spec1d_flux[spec1d.spec1d_mask==False],color='y',lw=0.3)

        dev = y-continuum0(spec1d.spec1d_pixel)
        outlier = np.where(dev>cr_rejection_tol*rms0)[0]
    
        y.mask[outlier] = True
        y.mask[outlier-cr_rejection_collateral] = True
        y.mask[outlier+cr_rejection_collateral] = True
#    plt.plot(spec1d.spec1d_pixel,continuum0(spec1d.spec1d_pixel),color='r')
#    plt.scatter(spec1d.spec1d_pixel[y.mask==True],spec1d.spec1d_flux[y.mask==True],color='magenta',s=20)
#    plt.plot(spec1d.spec1d_pixel[y.mask==False],spec1d.spec1d_flux[y.mask==False],color='k',lw=0.2)
#    print((spec1d.spec1d_pixel[y.mask==True],spec1d.spec1d_flux[y.mask==True]))
#    plt.show()
#    plt.close()
    return extract1d(aperture=spec1d.aperture,spec1d_pixel=spec1d.spec1d_pixel,spec1d_flux=spec1d.spec1d_flux,
                     spec1d_uncertainty=spec1d.spec1d_uncertainty,spec1d_mask=y.mask)


def column_stack(data,col):
    import numpy as np
    import mycode
    from ccdproc import Combiner
    import astropy.units as u
    import specutils
    from astropy.modeling import models
    from specutils.fitting import fit_lines
    from specutils.spectra import Spectrum1D
    import matplotlib.pyplot as plt

    column = data.data[:,col]
    column_uncertainty = data.uncertainty[:,col]
    column_mask = data.mask[:,col]

    stack = []
    sig_stack = []
    for j in range(0,len(column[0])):
        stack.append(data[:,col[j]])
        sig_stack.append(data.uncertainty._array[:,col[j]])
    sig_stack = np.array(sig_stack)
    c = Combiner(stack)
    c.weights = 1./sig_stack**2
    comb = c.average_combine(uncertainty_func=mycode.stdmean)
    lamb = np.arange(len(comb.data),dtype='float')  # unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    spec1d = Spectrum1D(flux=np.array(comb.data)*u.electron,spectral_axis=lamb*u.AA,uncertainty=comb.uncertainty,mask=comb.mask)
    return spec1d


def plot_apertures(pixel,spec,continuum,subregion,fit,realvirtual):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from specutils.manipulation import extract_region
    from astropy import units as u

    gs = plt.GridSpec(10,10)   #  define multi-panel plot
    gs.update(wspace=0,hspace=0)   #  specify inter-panel spacing
#    fig = plt.figure(figsize=(6,6)) # define plot size
    ax1 = fig.add_subplot(111)

    x = np.arange(0,len(spec))
    ax1.plot(pixel,spec,lw=0.3,color='k')
    ax1.plot(pixel,continuum(x),color='green',lw=1)
    ax1.set_xlim([0,np.max(x)])

    for j in range(0,len(fit)):
        use = np.where((pixel.value>=subregion[j].lower.value)&(pixel.value<=subregion[j].upper.value))[0]
        sub_spectrum_pixel = pixel[use]
        y_fit = fit[j](sub_spectrum_pixel).value+continuum(sub_spectrum_pixel.value)

        ax1.axvline(x=fit[j].mean.value,color='b',lw=0.3,alpha=0.7)
        ax1.plot(sub_spectrum_pixel.value,y_fit,color='r',lw=0.2,alpha=0.7)
        ax1.text(fit[j].mean.value,0,str(j+1),fontsize=8)
        if not realvirtual[j]:
            ax1.axvline(x=fit[j].mean.value,color='k',lw=0.5,alpha=1,linestyle='--')
    return ax1


def fit_aperture(spec,window,x_center):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from specutils.manipulation import extract_region
    from specutils import SpectralRegion
    from specutils.manipulation import extract_region
    import astropy.units as u
    from specutils.spectra import Spectrum1D
    from specutils.fitting import estimate_line_parameters
    from astropy.modeling import models,fitting
    from specutils.fitting import fit_lines

    val1 = x_center-window/2.  # window for fitting aperture
    val2 = x_center+window/2.  # window for fitting aperture
    subregion0 = SpectralRegion(val1*u.AA,val2*u.AA)  # define extraction region from window
    sub_spectrum0 = extract_region(spec,subregion0)  # extract from window
#    sub_spectrum = extract_region(spec,subregion0)#extract from window

    sub_spectrum = Spectrum1D(flux=sub_spectrum0.flux,spectral_axis=sub_spectrum0.spectral_axis,
                              uncertainty=sub_spectrum0.uncertainty)
#    print(sub_spectrum)
#    print(sub_spectrum)
    rough = estimate_line_parameters(sub_spectrum,models.Gaussian1D())  # get rough estimate of gaussian parameters for aperture
#    np.random.seed(0)
#    x = np.linspace(0.,10.,2056)
#    y = 3*np.exp(-0.5*(x-6.3)**2/0.8**2)
#    y+=np.random.normal(0.,0.2,x.shape)
#    spectrum = Spectrum1D(flux=y*u.Jy,spectral_axis=x*u.um)
#    g_init = models.Gaussian1D(amplitude=3.*u.Jy,mean=6.1*u.um,stddev=1.*u.um)
#    g_fit = fit_lines(spectrum,g_init)
#    y_fit = g_fit(x*u.um)
#    print(spectrum)

    g_init = models.Gaussian1D(amplitude=rough.amplitude.value*u.electron,
                               mean=rough.mean.value*u.AA,stddev=rough.stddev.value*u.AA)  # now do a fit using rough estimate as first guess
    g_fit0 = fit_lines(sub_spectrum,g_init)  # now do a fit using rough estimate as first guess
    return subregion0,g_fit0


def aperture_order(subregion,g_fit,realvirtual,initial):  # sort apertures in order of increasing column center
    import numpy as np

    means = [g_fit[i].mean.value for i in range(0,len(g_fit))]
    order0 = np.argsort(means)
    g_fit2 = [g_fit[order0[i]] for i in range(0,len(g_fit))]
    g_fit = g_fit2
    subregion2 = [subregion[order0[i]] for i in range(0,len(g_fit))]
    subregion = subregion2
    realvirtual2 = [realvirtual[order0[i]] for i in range(0,len(g_fit))]
    realvirtual = realvirtual2
    initial2 = [initial[order0[i]] for i in range(0,len(g_fit))]
    initial = initial2
    return subregion,g_fit,realvirtual,initial


def get_trace_xyz(j,columnspec_array,apertures_profile_middle,middle_column,trace_shift_max,trace_nlost_max):
    import numpy as np
    mean0 = apertures_profile_middle.fit[j].mean.value
    trace_x = []
    trace_y = []
    trace_z = []
    nlost = 0
    for k in reversed(range(0,middle_column)):
        if nlost<=trace_nlost_max:
            means = np.array([columnspec_array[k].apertures_profile.fit[q].mean.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
            keep = np.where(means==means)[0]
            means = means[keep]  # remove any NaNs due to poor fits
            amplitudes = np.array([columnspec_array[k].apertures_profile.fit[q].amplitude.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
            amplitudes = amplitudes[keep]
            dist = (mean0-means)**2
            if len(keep)>0:
                best = np.where(dist==np.min(dist))[0][0]
                if dist[best]<trace_shift_max:
                    nlost = 0
                    mean0 = means[best]
                    amplitude0 = amplitudes[best]
                    trace_x.append(np.median(columnspec_array[k].columns))
                    trace_y.append(means[best])
                    trace_z.append(amplitudes[best])
                else:
                    nlost = nlost+1
            else:
                nlost = nlost+1
    nlost = 0
    mean0 = apertures_profile_middle.fit[j].mean.value
    for k in range(middle_column,len(columnspec_array)):
        if nlost<trace_nlost_max:
            means = np.array([columnspec_array[k].apertures_profile.fit[q].mean.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
            keep = np.where(means==means)[0]
            means = means[keep]  # remove any NaNs due to poor fits
            amplitudes = np.array([columnspec_array[k].apertures_profile.fit[q].amplitude.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
            amplitudes = amplitudes[keep]
            dist = (mean0-means)**2
            if len(keep)>0:
                best = np.where(dist==np.min(dist))[0][0]
                if dist[best]<trace_shift_max:
                    nlost = 0
                    mean0 = means[best]
                    amplitude0 = amplitudes[best]
                    trace_x.append(np.median(columnspec_array[k].columns))
                    trace_y.append(means[best])
                    trace_z.append(amplitudes[best])
                else:
                    nlost = nlost+1
            else:
                nlost = nlost+1
    trace_x = np.array(trace_x)
    trace_y = np.array(trace_y)
    trace_z = np.array(trace_z)
    order = np.argsort(trace_x)
    trace_x = trace_x[order]
    trace_y = trace_y[order]
    trace_z = trace_z[order]
    return trace_x,trace_y,trace_z


def get_aperture(j,columnspec_array,apertures_profile_middle,middle_column,trace_order,trace_rejection_sigma,
                 trace_rejection_iterations,image_boundary,trace_shift_max,trace_nlost_max,
                 profile_rejection_iterations,profile_nsample,profile_order,window):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models,fitting
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D

    trace_x,trace_y,trace_z = get_trace_xyz(j,columnspec_array,apertures_profile_middle,
                                            middle_column,trace_shift_max,trace_nlost_max)
    trace_func_init = models.Polynomial1D(degree=trace_order)
    fitter = fitting.LinearLSQFitter()
    y = np.ma.masked_array(trace_y,mask=False)
    z = np.ma.masked_array(trace_z,mask=False)
    sig_y = trace_y  # fudge for now, not used for anything

    trace_pix_min = []
    trace_pix_max = []
    for i in range(0,len(trace_x)):
        trace_pix_min.append(np.interp(x=trace_y[i],xp=[image_boundary.lower[q][1] for q in range(0,len(image_boundary.lower))],fp=[image_boundary.lower[q][0] for q in range(0,len(image_boundary.lower))]))
        trace_pix_max.append(np.interp(x=trace_y[i],xp=[image_boundary.upper[q][1] for q in range(0,len(image_boundary.upper))],fp=[image_boundary.upper[q][0] for q in range(0,len(image_boundary.upper))]))
    trace_pix_min = np.array(trace_pix_min)
    trace_pix_max = np.array(trace_pix_max)

    y.mask[trace_x<trace_pix_min] = True
    y.mask[trace_x>trace_pix_max] = True
    z.mask[trace_x<trace_pix_min] = True
    z.mask[trace_x>trace_pix_max] = True

    trace_func = trace_func_init
    trace_rms = np.sqrt(np.mean((y-trace_func(trace_x))**2))

    for q in range(0,trace_rejection_iterations):
        use = np.where(y.mask==False)[0]
        if len(use)>=3:
            trace_func = fitter(trace_func_init,trace_x[use],trace_y[use])
            trace_rms = np.sqrt(np.mean((y[use]-trace_func(trace_x)[use])**2))
            outlier = np.where(np.abs(trace_y-trace_func(trace_x))>trace_rejection_sigma*trace_rms)[0]
            y.mask[outlier] = True
    use = np.where(y.mask==False)[0]
    trace_npoints = len(use)

    order = [trace_order]
    rejection_sigma = [trace_rejection_sigma]
#    pix_min = [trace_pix_min]
#    pix_max = [trace_pix_max]
    func = [trace_func]
    rms = [trace_rms]
    npoints = [trace_npoints]
    if len(trace_x[y.mask==False])>0:
        pix_min = [np.min(trace_x[y.mask==False])]
        pix_max = [np.max(trace_x[y.mask==False])]
    else:
        pix_min = [np.float(-999.)]
        pix_max = [np.float(-999.)]


    print(j,trace_npoints,pix_min,pix_max)

#    fig=plt.figure(1)
#    ax1,ax2=plot_trace(j,trace_x,trace_y,trace_z,y,trace_func,np.median(trace_pix_min),np.median(trace_pix_max),trace_rms,trace_npoints,trace_order,fig)

#    print('press \'d\' to delete point nearest cursor')
#    print('press \'c\' to un-delete deleted point nearest cursor')
#    print('press \'l\' to re-set lower limit of fit region')
#    print('press \'u\' to re-set lower limit of fit region')
#    print('press \'i\' to iterate outlier rejection')
#    print('press \'r\' to change rejection sigma threshold')
#    print('press \'o\' to change order of fit')
#    print('press \'z\' to return to original and erase edits')
#    cid=fig.canvas.mpl_connect('key_press_event',lambda event: on_key_trace(event,[j,trace_x,trace_y,trace_z,y,sig_y,pix_min,pix_max,func,rms,npoints,order,rejection_sigma,trace_rejection_iterations,trace_func_init,fig]))
#    plt.show()
#    plt.close()

    if len(trace_x[y.mask==False])>0:
        pix_min.append(np.min(trace_x[y.mask==False]))
        pix_max.append(np.max(trace_x[y.mask==False]))
    else:
        pix_min.append(np.float(-999.))
        pix_max.append(np.float(-999.))
    profile_init = models.Polynomial1D(degree=profile_order)
    profile_fitter = fitting.LinearLSQFitter()

    x = np.array([np.median(columnspec_array[k].columns) for k in range(0,len(columnspec_array))])
    if apertures_profile_middle.realvirtual[j]:  # if this is a real aperture, fit line profile as function of pixel
        y1 = []
        y2 = []
        for i in range(0,len(x)):  # re-fit Gaussian profile at different x values along spectral axis
            center = func[len(func)-1](x[i])
            if((center>0.)&(center<np.max(columnspec_array[i].pixel.value))):  # make sure the trace function y(x) makes sense at this x
                spec1d = Spectrum1D(spectral_axis=columnspec_array[i].pixel,flux=columnspec_array[i].spec*u.electron,
                                    uncertainty=columnspec_array[i].err,mask=columnspec_array[i].mask)
                subregion0,g_fit0 = fit_aperture(spec1d-columnspec_array[i].continuum(columnspec_array[i].pixel.value),window,center)
                y1.append(g_fit0.stddev.value)
                y2.append(g_fit0.amplitude.value)
            else:  # otherwise give a place-holder value and mask it below
                y1.append(np.float(-999.))
                y2.append(np.float(-999.))
        y1 = np.array(y1)
        y2 = np.array(y2)
        sigma_y = np.ma.masked_array(y1,mask=False)
        sigma_y.mask[y1<-998.] = True
        sigma_y.mask[x<pix_min[len(pix_min)-1]] = True
        sigma_y.mask[x>pix_max[len(pix_max)-1]] = True
        amplitude_y = np.ma.masked_array(y2,mask=False)
        amplitude_y.mask[y2<-998.] = True
        amplitude_y.mask[x<pix_min[len(pix_min)-1]] = True
        amplitude_y.mask[x>pix_max[len(pix_max)-1]] = True

        profile_sigma = profile_init
        profile_sigma_rms = np.sqrt(np.mean((y[use]-profile_sigma(x[use]))**2))
        profile_amplitude = profile_init
        profile_amplitude_rms = np.sqrt(np.mean((y[use]-profile_amplitude(x[use]))**2))
        for q in range(0,profile_rejection_iterations):
            use = np.where((sigma_y.mask==False)&(amplitude_y.mask==False))[0]
            if len(use)>3.:
                profile_sigma = profile_fitter(profile_init,x[use],sigma_y[use])
                profile_sigma_rms = np.sqrt(np.mean((sigma_y[use]-profile_sigma(x[use]))**2))
                profile_amplitude = profile_fitter(profile_init,x[use],amplitude_y[use])
                profile_amplitude_rms = np.sqrt(np.mean((amplitude_y[use]-profile_amplitude(x[use]))**2))
                outlier1 = (np.where(np.abs(sigma_y-profile_sigma(x))>3.*profile_sigma_rms))[0]  # reject outliers
                outlier2 = (np.where(np.abs(amplitude_y-profile_amplitude(x))>3.*profile_amplitude_rms))[0]  # reject outliers
                sigma_y.mask[outlier1] = True
                amplitude_y.mask[outlier2] = True
        use = np.where((sigma_y.mask==False)&(amplitude_y.mask==False))[0]
        profile_npoints = len(use)

#        print(profile_sigma_rms/np.median(sigma_y[sigma_y.mask==False]),profile_amplitude_rms/np.median(amplitude_y[amplitude_y.mask==False]))
#        plt.scatter(x,amplitude_y,s=3,color='k')
#        plt.scatter(x[amplitude_y.mask==True],amplitude_y[amplitude_y.mask==True],s=13,marker='x',color='y')
#        plt.plot(x,profile_amplitude(x),color='r')
#        plt.axvline(x=pix_min[len(pix_min)-1],lw=0.3,linestyle=':',color='k')
#        plt.axvline(x=pix_max[len(pix_max)-1],lw=0.3,linestyle=':',color='k')
##        plt.ylim([1,3])
#        plt.xlabel('x [pixel]')
#        plt.ylabel('LSF sigma [pixel]')
#        plt.show()
#        plt.close()
    else:
        profile_sigma = -999
        profile_sigma_rms = -999
        profile_amplitude = -999
        profile_amplitude_rms = -999
        profile_npoints = -999
    aperture0 = aperture(trace_aperture=j+1,trace_func=func[len(func)-1],trace_rms=rms[len(rms)-1],
                         trace_npoints=npoints[len(npoints)-1],trace_pixel_min=pix_min[len(pix_min)-1],
                         trace_pixel_max=pix_max[len(pix_max)-1],profile_sigma=profile_sigma,
                         profile_sigma_rms=profile_sigma_rms,profile_amplitude=profile_amplitude,
                         profile_amplitude_rms=profile_amplitude_rms,profile_npoints=profile_npoints)
    return aperture0


class mapflat_compare:
    def __init__(self,func=None,rms=None,npoints=None):
        self.func=func
        self.rms=rms
        self.npoints=npoints

class wavcal:
    def __init__(self,wav=None,nthar=None,thar_npoints=None,thar_rms=None,wav_rms=None,mjd=None,line_sigma=None):
        self.wav = wav
        self.nthar = nthar
        self.thar_npoints=thar_npoints
        self.thar_rms = thar_rms
        self.wav_rms = wav_rms
        self.mjd = mjd
        self.line_sigma = line_sigma

class scatteredlightfunc:
    def __init__(self,func=None,rms=None):
        self.func = func
        self.rms = rms

class image_boundary:
    def __init__(self,lower=None,upper=None):
        self.lower = lower
        self.upper = upper

class linelist:
    def __init__(self,wavelength=None,species=None):
        self.wavelength = wavelength
        self.species = species

class aperture_profile:
    def __init__(self,fit=None,subregion=None,realvirtual=None,initial=None):
        self.fit = fit
        self.subregion = subregion
        self.realvirtual = realvirtual
        self.initial = initial

class aperture:
    def __init__(self,trace_aperture=None,trace_func=None,trace_rms=None,trace_npoints=None,trace_pixel_min=None,trace_pixel_max=None,
                 profile_sigma=None,profile_sigma_rms=None,profile_amplitude=None,profile_amplitude_rms=None,profile_npoints=None):
        self.trace_aperture = trace_aperture
        self.trace_func = trace_func
        self.trace_rms = trace_rms
        self.trace_npoints = trace_npoints
        self.trace_pixel_min = trace_pixel_min
        self.trace_pixel_max = trace_pixel_max
        self.profile_sigma = profile_sigma
        self.profile_sigma_rms = profile_sigma_rms
        self.profile_amplitude = profile_amplitude
        self.profile_amplitude_rms = profile_amplitude_rms
        self.profile_npoints = profile_npoints

class thars:
    def __init__(self,aperture=None,wav=None,pix=None,wav_func=None,vel_func=None,pix_std=None,wav_func_std=None,vel_func_std=None,temperature=None):
        self.aperture = aperture
        self.wav = wav
        self.pix = pix
        self.wav_func = wav_func
        self.vel_func = vel_func
        self.pix_std = pix_std
        self.wav_func_std = wav_func_std
        self.vel_func_std = vel_func_std
        self.temperature = temperature


def get_aperture_profile(apertures_initial,spec1d,continuum,window):
    import numpy as np

    subregion = []
    g_fit = []
    realvirtual = []  # 1 for real aperture, 0 for virtual aperture (ie a placeholder aperture for bad/unplugged fiber etc.)
    initial = []  # 1 for apertures identified automatically, 0 for those identified by hand

    for j in range(0,len(apertures_initial)):
        x_center = apertures_initial['line_center'][j].value
        subregion0,g_fit0 = fit_aperture(spec1d-continuum(spec1d.spectral_axis.value),window,x_center)
        subregion.append(subregion0)
        g_fit.append(g_fit0)
        realvirtual.append(True)
        initial.append(True)
    return aperture_profile(g_fit,subregion,realvirtual,initial)


class thar:
    def __init__(self,aperture=None,wav=None,pix=None,hjd=None):
        self.aperture = aperture
        self.wav = wav
        self.pix = pix
        self.hjd = hjd

class columnspec:
    def __init__(self,columns=None,spec=None,mask=None,err=None,pixel=None,continuum=None,rms=None,
                 apertures_initial=None,apertures_profile=None):
        self.columns = columns
        self.spec = spec
        self.mask = mask
        self.err = err
        self.pixel = pixel
        self.continuum = continuum
        self.rms = rms
        self.apertures_initial = apertures_initial
        self.apertures_profile = apertures_profile

class extract1d:
    def __init__(self,aperture=None,spec1d_pixel=None,spec1d_flux=None,spec1d_uncertainty=None,spec1d_mask=None):
        self.aperture = aperture
        self.spec1d_pixel = spec1d_pixel
        self.spec1d_flux = spec1d_flux
        self.spec1d_uncertainty = spec1d_uncertainty
        self.spec1d_mask = spec1d_mask

class id_lines:
    def __init__(self,aperture=None,fit_lines=None,pix=None,wav=None,func=None,rms=None,npoints=None,
                 resolution=None,resolution_rms=None,resolution_npoints=None):
        self.aperture = aperture
        self.fit_lines = fit_lines
        self.wav = wav
        self.func = func
        self.rms = rms
        self.npoints = npoints
        self.resolution = resolution
        self.resolution_rms = resolution_rms
        self.resolution_npoints = resolution_npoints

def get_columnspec(data,trace_step,n_lines,continuum_rejection_low,continuum_rejection_high,
                   continuum_rejection_iterations,continuum_rejection_order,threshold_factor,window):
    import numpy as np
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative

    n_cols = np.shape(data)[1]
    trace_n = np.long(n_cols/trace_step)
#    print(n_cols,trace_n)
    trace_cols = np.linspace(0,n_cols,trace_n,dtype='int')

#    apertures_initial = []
#    col = []
#    spec1d = []
#    continuum = []
#    rms = []
#    pixel = []
#    apertures_profile = []
    columnspec_array = []
    for i in range(0,len(trace_cols)-1):
        print('working on '+str(i+1)+' of '+str(len(trace_cols))+' trace columns')
        col0 = np.arange(n_lines)+trace_cols[i]
        spec1d0 = column_stack(data,col0)
#        np.pause()
        continuum0,rms0 = get_continuum(spec1d0,continuum_rejection_low,continuum_rejection_high,
                                        continuum_rejection_iterations,continuum_rejection_order)
        pixel0 = (np.arange(len(spec1d0.data),dtype='float'))*u.AA  # unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
        spec_contsub = spec1d0-continuum0(pixel0.value)
        spec_contsub.uncertainty.quantity.value[:] = rms0
#        spec_contsub.mask[:] = False
        apertures_initial0 = find_lines_derivative(spec_contsub,flux_threshold=threshold_factor*rms0)  # find peaks in continuum-subtracted "spectrum"
#        apertures_profile0 = get_aperture_profile(spec1d0,window)

        apertures_profile0 = get_aperture_profile(apertures_initial0,spec1d0,continuum0,window)
        columnspec0 = columnspec(columns=col0,spec=spec1d0.data,mask=spec1d0.mask,
                                 err=spec1d0.uncertainty,pixel=pixel0,continuum=continuum0, # have to break up spec1d into individual
                                 rms=rms0,apertures_initial=apertures_initial0,             # data/mask/uncertainty fields in order
                                 apertures_profile=apertures_profile0)                      # to allow pickling of columnspec_array
        print('found '+str(len(apertures_initial0))+' apertures in column '+str(col0))
        columnspec_array.append(columnspec0)

    return columnspec_array


def fiddle_apertures(columnspec_array,column,window,apertures,find_apertures_file):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models

    subregion,fit,realvirtual,initial = apertures.subregion,apertures.fit,apertures.realvirtual,apertures.initial
    subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,apertures.initial)

    fig = plt.figure(1)
    x = np.arange(0,len(columnspec_array[column].spec))
    ax1 = fig.add_subplot(111)
    ax1.plot(columnspec_array[column].pixel,columnspec_array[column].spec,lw=0.3,color='k')
    ax1.plot(columnspec_array[column].pixel,columnspec_array[column].continuum(x),color='green',lw=1)
    ax1.set_xlim([0,np.max(x)])
    for j in range(0,len(fit)):
        use = np.where((columnspec_array[column].pixel.value>=subregion[j].lower.value)&(columnspec_array[column].pixel.value<=subregion[j].upper.value))[0]
        sub_spectrum_pixel = columnspec_array[column].pixel[use]
        y_fit = fit[j](sub_spectrum_pixel).value+columnspec_array[column].continuum(sub_spectrum_pixel.value)

        ax1.axvline(x=fit[j].mean.value,color='b',lw=0.3,alpha=0.7)
        ax1.plot(sub_spectrum_pixel.value,y_fit,color='r',lw=0.2,alpha=0.7)
        ax1.set_ylim([0,10000])
        ax1.text(fit[j].mean.value,0,str(j+1),fontsize=8)
        if not realvirtual[j]:
            ax1.axvline(x=fit[j].mean.value,color='k',lw=0.5,alpha=1,linestyle='--')

    print('press \'d\' to delete aperture nearest cursor \n')
    print('press \'e\' to delete all apertures \n')
    print('press \'n\' to add new real aperture at cursor position \n')
    print('press \'a\' to add new phantom aperture at cursor position \n')
    print('press \'z\' to return to initial apertures \n')
    print('press \'q\' to quit \n')
    cid = fig.canvas.mpl_connect('key_press_event',lambda event: on_key_find(event,[columnspec_array,column,
                                                                                    subregion,fit,realvirtual,initial,window,fig]))
    plt.show()
    subregion,fit,realvirtual,initial = aperture_order(subregion,fit,realvirtual,initial)
    return aperture_profile(fit,subregion,realvirtual,initial)


def get_image_boundary(data,image_boundary_fiddle,image_boundary0):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from specutils import SpectralRegion
    import astropy.units as u
    from astropy.modeling import models

    lower_x = []
    lower_y = []
    upper_x = []
    upper_y = []
    reject_x = []
    reject_y = []
    if image_boundary_fiddle:
        for i in range(0,len(image_boundary0.lower)):
            x,y = image_boundary0.lower[i]
            lower_x.append(x)
            lower_y.append(y)
        for i in range(0,len(image_boundary0.upper)):
            x,y = image_boundary0.upper[i]
            upper_x.append(x)
            upper_y.append(y)
#        for i in range(0,len(image_boundary0.reject)):
#            x,y = image_boundary0.reject[i]
#            reject_x.append(x)
#            reject_y.append(y)

    fig = plt.figure(1)
    ax1 = plot_boundary(lower_x,lower_y,upper_x,upper_y,reject_x,reject_y,data,fig)
    print('click \'l\' to add a lower boundary point')
    print('click \'u\' to add an upper boundary point')
    print('click \'b\' to mark rejection line')
    print('click \'q\' when finished')
    cid = fig.canvas.mpl_connect('key_press_event',lambda event: on_key_boundary(event,[data,lower_x,lower_y,upper_x,
                                                                                        upper_y,reject_x,reject_y,fig]))
    plt.show()
    plt.close()

    x1 = np.array(lower_x)
    y1 = np.array(lower_y)
    order = np.argsort(y1)
    x1 = x1[order]
    y1 = y1[order]
    x1 = np.append(x1[0],x1)
    y1 = np.append(0,y1)
    x1 = np.append(x1,x1[len(x1)-1])
    y1 = np.append(y1,len(data.data))

    x2 = np.array(upper_x)
    y2 = np.array(upper_y)
    order = np.argsort(y2)
    x2 = x2[order]
    y2 = y2[order]
    x2 = np.append(x2[0],x2)
    y2 = np.append(0,y2)
    x2 = np.append(x2,x2[len(x2)-1])
    y2 = np.append(y2,len(data.data))
    
    return image_boundary(lower=[(x1[q],y1[q]) for q in range(0,len(x1))],upper=[(x2[q],y2[q]) for q in range(0,len(x2))])


def get_id_lines_template(extract1d,linelist,continuum_rejection_low,continuum_rejection_high,
                          continuum_rejection_iterations,continuum_rejection_order,threshold_factor,
                          window,id_lines_order,id_lines_rejection_iterations,id_lines_rejection_sigma,
                          id_lines_tol_angs,id_lines_template_fiddle,id_lines_template0,resolution_order,
                          resolution_rejection_iterations):
    import numpy as np
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative
    from specutils.spectra import Spectrum1D
    import matplotlib.pyplot as plt
    from copy import deepcopy
    from astropy.modeling import models,fitting

    continuum0,spec_contsub,fit_lines = get_fitlines(extract1d,continuum_rejection_low,continuum_rejection_high,
                                                     continuum_rejection_iterations,threshold_factor,window,
                                                     continuum_rejection_order)

    line_centers = [fit_lines.fit[j].mean.value for j in range(0,len(fit_lines.fit))]
    id_lines_pix = []
    id_lines_wav = []
    id_lines_species = []
    id_lines_used = []
    order = [id_lines_order]
    rejection_iterations = [id_lines_rejection_iterations]
    rejection_sigma = [id_lines_rejection_sigma]
    func = [models.Legendre1D(degree=1)]
    rms = []
    npoints = []
    if id_lines_template_fiddle:
        id_lines_used = np.where(id_lines_template0.wav.mask==False)[0].tolist()
        id_lines_pix = [id_lines_template0.fit_lines.fit[q].mean.value for q in id_lines_used]
        id_lines_wav = id_lines_template0.wav[id_lines_used].tolist()
        func = [id_lines_template0.func]
        order = [id_lines_template0.func.degree]
        rms = [id_lines_template0.rms]
        npoints = [id_lines_template0.npoints]

    fig = plt.figure(1)
    ax1,ax2 = plot_id_lines(extract1d,continuum0,fit_lines,line_centers,id_lines_pix,id_lines_wav,
                            id_lines_used,func[len(func)-1],fig)
#
    print('press \'m\' to ID line nearest cursor \n')
    print('press \'d\' to delete ID for line nearest cursor \n')
    print('press \'o\' to change order of polynomial \n')
    print('press \'r\' to change rejection sigma factor \n')
    print('press \'t\' to change number of rejection iterations \n')
    print('press \'g\' to re-fit wavelength solution \n')
    print('press \'l\' to add lines from linelist according to fit \n')
    print('press \'q\' to quit \n')
    print('press \'.\' to print cursor position and position of nearest line \n')

    cid = fig.canvas.mpl_connect('key_press_event',lambda event: on_key_id_lines(event,[deepcopy(extract1d),continuum0,
                                                                                        fit_lines,linelist,line_centers,
                                                                                        id_lines_pix,id_lines_wav,id_lines_used,
                                                                                        order,rejection_iterations,rejection_sigma,
                                                                                        func,rms,npoints,id_lines_tol_angs,fig]))
    plt.show()
    plt.close()

    wav = np.ma.masked_array(np.full((len(fit_lines.fit)),-999,dtype='float'),
                             mask=np.full((len(fit_lines.fit)),True,dtype=bool))

    for j in range(0,len(id_lines_pix)):
        wav[id_lines_used[j]] = id_lines_wav[j]
        wav.mask[id_lines_used[j]] = False

    resolution,resolution_rms,resolution_npoints = get_resolution(deepcopy(fit_lines),deepcopy(wav),
                                                                  resolution_order,resolution_rejection_iterations)

    return id_lines(aperture=extract1d.aperture,fit_lines=fit_lines,wav=wav,func=func[len(func)-1],
                    rms=rms[len(rms)-1],npoints=npoints[len(npoints)-1],resolution=resolution,resolution_rms=resolution_rms,
                    resolution_npoints=resolution_npoints)


def get_fitlines(extract1d,continuum_rejection_low,continuum_rejection_high,continuum_rejection_iterations,
                 threshold_factor,window,continuum_rejection_order):
    import numpy as np
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative
    from copy import deepcopy
    import matplotlib.pyplot as plt

    spec1d0 = Spectrum1D(spectral_axis=deepcopy(extract1d.spec1d_pixel)*u.AA,flux=deepcopy(extract1d.spec1d_flux),
                         uncertainty=deepcopy(extract1d.spec1d_uncertainty),mask=deepcopy(extract1d.spec1d_mask))
    continuum0,rms0 = get_continuum(spec1d0,continuum_rejection_low,continuum_rejection_high,
                                    continuum_rejection_iterations,continuum_rejection_order)
    pixel0 = extract1d.spec1d_pixel*u.AA
    spec_contsub = spec1d0-continuum0(pixel0.value)
    spec_contsub.uncertainty.quantity.value[:] = rms0
    spec_contsub.mask = spec1d0.mask
    spec_contsub.flux[spec_contsub.flux!=spec_contsub.flux] = 0.
#    spec_contsub.flux[spec_contsub.flux.value<-100] = 0.*u.electron
    id_lines_initial0 = find_lines_derivative(spec_contsub,flux_threshold=threshold_factor*rms0)  # find peaks in continuum-subtracted "spectrum"
    fit_lines = get_aperture_profile(id_lines_initial0,spec1d0,continuum0,window)
    return continuum0,spec_contsub,fit_lines


def get_id_lines_translate(extract1d_template,id_lines_template,extract1d,linelist,continuum_rejection_low,
                           continuum_rejection_high,continuum_rejection_iterations,continuum_rejection_order,
                           threshold_factor,window,id_lines_order,id_lines_rejection_iterations,
                           id_lines_rejection_sigma,id_lines_tol_angs,id_lines_tol_pix,resolution_order,
                           resolution_rejection_iterations,add_lines_iterations):
    import numpy as np
    import scipy
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative
    from specutils.spectra import Spectrum1D
    import matplotlib.pyplot as plt
    from copy import deepcopy
    from astropy.modeling import models,fitting
    from dynesty import DynamicNestedSampler,NestedSampler

    fit_lines = [0]
    wav = [0.]
    func = [0]
    rms = [-999.]
    npoints = [0]
    resolution = 0.
    resolution_rms = -999.
    resolution_npoints = 0

    print(len(np.where(extract1d.spec1d_mask==False)[0]),'asfdasdfasdfasdf')    
    if len(np.where(extract1d.spec1d_mask==False)[0])>100:
#        plt.plot(extract1d.spec1d_flux[extract1d.spec1d_mask==False])
#        plt.show()
#        plt.close()

        continuum0,spec_contsub,fit_lines = get_fitlines(extract1d,continuum_rejection_low,continuum_rejection_high,
                                                         continuum_rejection_iterations,threshold_factor,window,
                                                         continuum_rejection_order)
        pixelmin = np.min(extract1d.spec1d_pixel[extract1d.spec1d_mask==False])
        pixelmax = np.max(extract1d.spec1d_pixel[extract1d.spec1d_mask==False])
        pixelscale = 0.5*(pixelmax-pixelmin)
        pixel0 = pixelmin+pixelscale
        use = np.where((extract1d.spec1d_pixel>=pixelmin)&(extract1d.spec1d_pixel<=pixelmax)&(extract1d.spec1d_mask==False)&(extract1d.spec1d_flux.value<1.e+6))[0]
#        use=np.where((extract1d.spec1d_pixel>=pixelmin)&(extract1d.spec1d_pixel<=pixelmax)&(extract1d.spec1d_mask==False))[0]

        continuum0_template,spec_contsub_template,fit_lines_template = get_fitlines(extract1d_template,continuum_rejection_low,
                                                                                    continuum_rejection_high,continuum_rejection_iterations,
                                                                                    threshold_factor,window,continuum_rejection_order)
        pixelmin_template = np.min(extract1d_template.spec1d_pixel[extract1d_template.spec1d_mask==False])
        pixelmax_template = np.max(extract1d_template.spec1d_pixel[extract1d_template.spec1d_mask==False])
        pixelscale_template = 0.5*(pixelmax_template-pixelmin_template)
        pixel0_template = pixelmin_template+pixelscale_template
        use_template = np.where((extract1d_template.spec1d_pixel>=pixelmin_template)&(extract1d_template.spec1d_pixel<=pixelmax_template)&(extract1d.spec1d_mask==False))[0]

        if len(fit_lines.fit)>0:

#        print(np.median(spec_contsub0.flux.value[use]))
#        spec_contsub=spec_contsub0/(np.median(spec_contsub0.flux.value[use]))
#        spec_contsub_template=spec_contsub_template0/(np.median(spec_contsub_template0.flux.value[use_template]))
            '''first run dynesty to fit for zero order offset between current lamp spectrum and template (gradient-descent gets stuck in local minima)'''
            ndim = 8
            def get_interp(q):  # spec_contsub,spec_contsub_template,use,use_template,pixel0_template,pixelscale_template):
                import numpy as np
                import matplotlib.pyplot as plt
        #    spec_contsub,spec_contsub_template,use,use_template,pixel0_template,pixelscale_template = args
                x0 = np.array(spec_contsub_template.spectral_axis.value[use_template],dtype='float')
                y0 = spec_contsub_template.flux.value[use_template]
                xscale = (x0-pixel0_template)/pixelscale_template
                x1 = q[1]*pixelscale_template+x0*(1.+np.polynomial.polynomial.polyval(xscale,q[2:]))
#            x1 = q[1]*pixelscale_template+x0*(1.+np.polynomial.chebyshev.chebval(xscale,q[2:]))
                interp = q[0]*np.interp(spec_contsub.spectral_axis.value[use],x1,y0)
                return interp

            def my_func(q):  # spec_contsub,spec_contsub_template,use,use_template,pixel0_template,pixelscale_template):
                import numpy as np
                import matplotlib.pyplot as plt
                interp = get_interp(q)  # ,spec_contsub,spec_contsub_template,use,use_template,pixel0_template,pixelscale_template)
#        diff2 = ((spec_contsub.flux.value[use]-interp)/spec_contsub.uncertainty.quantity.value[use])**2
                diff2 = (spec_contsub.flux.value[use]-interp)**2
                return np.sum(diff2)

            def loglike(q):
                return -0.5*my_func(q)  # ,spec_contsub,spec_contsub_template,use,use_template,pixel0_template,pixelscale_template)

            def ptform(u):
                prior = []
                prior.append([0.,5.])
                prior.append([-2.0,2.0])
                prior.append([-0.0,0.0])
                prior.append([-0.0,0.0])
                prior.append([-0.0,0.0])
                prior.append([-0.0,0.0])
                prior.append([-0.0,0.0])
                prior.append([-0.0,0.0])
                prior = np.array(prior)
                x = np.array(u)
                for i in range(0,len(x)):
                    x[i] = prior[i][0]+(prior[i][1]-prior[i][0])*u[i]
                return x

            dsampler = DynamicNestedSampler(loglike,ptform,ndim,bound='multi')
#        sampler = NestedSampler(loglike,ptform,ndim,bound='multi')
            dsampler.run_nested(maxcall=12000)
#        sampler.run_nested(dlogz=0.05)
            best = np.where(dsampler.results.logl==np.max(dsampler.results.logl))[0][0]
#        best=np.where(sampler.results.logl==np.max(sampler.results.logl))[0][0]

            '''now run gradient descent to find higher-order corrections to get best match'''    
#        params = np.append(dsampler.results.samples[best],np.array([0.,0.,0.,0.,0.]))
            params = np.append(dsampler.results.samples[best],np.array([0.,0.,0.,0.,0.]))
#        params = np.array([q0,q1,0.,0.,0.,0.,0.])

#        print('params: ',params)
#        print('log likelihood = ',loglike(params)/1.e9)
            shiftstretch = scipy.optimize.minimize(my_func,params,method='Powell')
            interp = get_interp(shiftstretch.x)
            print('log likelihood/1e9 = ',loglike(shiftstretch.x)/1.e9)

#            plt.plot(spec_contsub.spectral_axis[use],spec_contsub.flux[use],color='g',lw=0.3,alpha=0.9)
#            plt.plot(spec_contsub_template.spectral_axis[use_template],
#                     spec_contsub_template.flux[use_template],color='r',lw=0.3,alpha=0.9)
#            plt.plot(spec_contsub_template.spectral_axis[use],interp,color='cyan',lw=0.3,alpha=0.9)
#            plt.show()
#            plt.close()

            found = 0
            id_lines_pix = []
            id_lines_wav = []
            id_lines_species = []
            id_lines_used = []
            order = [id_lines_template.func.degree]
            rejection_iterations = [id_lines_rejection_iterations]
            rejection_sigma = [id_lines_rejection_sigma]

            x0 = np.array([id_lines_template.fit_lines.fit[q].mean.value for q in range(0,len(id_lines_template.fit_lines.fit))])

            xscale = (x0-pixel0_template)/pixelscale_template
            x1 = np.float(shiftstretch.x[1])*pixelscale_template+x0*(1.+np.polynomial.polynomial.polyval(xscale,shiftstretch.x[2:]))
#        x1 = np.float(shiftstretch.x[1])*pixelscale_template+x0*(1.+np.polynomial.chebyshev.chebval(xscale,shiftstretch.x[2:]))
            for i in range(0,len(x0)):
                if id_lines_template.wav.mask[i]==False:
                    dist = np.sqrt((x1[i]-np.array([fit_lines.fit[q].mean.value for q in range(0,len(fit_lines.fit))]))**2)
                    nan = np.where(dist!=dist)[0]
                    if len(nan)>0:
                        dist[nan] = np.inf  # artificially replace any nans with huge values
                    best = np.where(dist==np.min(dist))[0][0]
#                    print(best,dist[best],id_lines_tol_pix)
                    if dist[best]<id_lines_tol_pix:
                        id_lines_pix.append(fit_lines.fit[best].mean.value)
                        id_lines_wav.append(id_lines_template.wav[i])
                        id_lines_used.append(best)

            if len(id_lines_pix)>0:
                func = [models.Legendre1D(degree=1)]
                rms = []
                npoints = []
                func0,rms0,npoints0,y = id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,
                                                     order,rejection_iterations,rejection_sigma)
                for j in range(0,add_lines_iterations):
                    new_pix,new_wav,new_used = line_id_add_lines(linelist,[fit_lines.fit[q].mean.value for q in range(0,len(fit_lines.fit))],id_lines_used,func0,id_lines_tol_angs)
                    for i in range(0,len(new_pix)):
                        id_lines_pix.append(new_pix[i])
                        id_lines_wav.append(new_wav[i])
                        id_lines_used.append(new_used[i])
                    func0,rms0,npoints0,y = id_lines_fit(id_lines_pix,id_lines_wav,id_lines_used,
                                                         order,rejection_iterations,rejection_sigma)
                func.append(func0)
                rms.append(rms0)
                npoints.append(npoints0)

                wav = np.ma.masked_array(np.full((len(fit_lines.fit)),-999,dtype='float'),mask=np.full((len(fit_lines.fit)),True,dtype=bool))
                for j in range(0,len(id_lines_pix)):
                    wav[id_lines_used[j]] = id_lines_wav[j]
                    wav.mask[id_lines_used[j]] = False

                resolution,resolution_rms,resolution_npoints = get_resolution(deepcopy(fit_lines),deepcopy(wav),
                                                                              resolution_order,resolution_rejection_iterations)
        
    return id_lines(aperture=extract1d.aperture,fit_lines=fit_lines,wav=wav,func=func[len(func)-1],rms=rms[len(rms)-1],
                    npoints=npoints[len(npoints)-1],resolution=resolution,resolution_rms=resolution_rms,
                    resolution_npoints=resolution_npoints)


#######

    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    ax1.set_xlim([np.min(extract1d.spec1d_pixel[use]),np.max(extract1d.spec1d_pixel[use])])
    ax1.plot(extract1d.spec1d_pixel[use],extract1d.spec1d_flux[use],color='k',lw=0.5)
    ax1.plot(extract1d.spec1d_pixel[use],continuum0(extract1d.spec1d_pixel[use]),color='b',lw=0.3)
    ax1.plot(extract1d.spec1d_pixel[use],interp,color='orange',lw=0.3)
#    ax1.plot(x1,y0,color='orange',lw=0.3)
    ax1.set_xlabel('pixel')
    ax1.set_ylabel('counts')
    for j in range(0,len(fit_lines.fit)):
        use = np.where((extract1d.spec1d_pixel>=fit_lines.subregion[j].lower.value)&(extract1d.spec1d_pixel<=fit_lines.subregion[j].upper.value))[0]
        sub_spectrum_pixel = extract1d.spec1d_pixel[use]*u.AA
        y_fit = fit_lines.fit[j](sub_spectrum_pixel).value+continuum0(sub_spectrum_pixel.value)
        ax1.plot(sub_spectrum_pixel,y_fit,color='r',lw=0.3)
        #            x_center=id_lines_initial0['line_center'][j].value
        #            print(fit_lines.fit[j].mean.value)
        plt.axvline(x=fit_lines.fit[j].mean.value,linestyle=':',color='0.3',lw=0.3)

    ax1.plot(extract1d_template.spec1d_pixel[extract1d_template.spec1d_mask==False],
             extract1d_template.spec1d_flux[extract1d_template.spec1d_mask==False],color='cyan',lw=0.3)
    
    line_centers = [fit_lines.fit[j].mean for j in range(0,len(fit_lines.fit))]
    id_lines_pix = []
    id_lines_wav = []
    id_lines_species = []
    id_lines_used = []
    order = [id_lines_order]
    rejection_iterations = [id_lines_rejection_iterations]
    rejection_sigma = [id_lines_rejection_sigma]
#        func = [models.Polynomial1D(degree=1)]
    func = [models.Legendre1D(degree=1)]
    rms = []
    npoints = []
    plt.show()
    plt.close()


def get_extract1d(j,data,apertures_profile_middle,aperture_array,aperture_peak,pix,extract1d_aperture_width,offset):
    import numpy as np
    import scipy
    from astropy.nddata import CCDData
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from astropy.nddata import StdDevUncertainty

    above0,below0 = get_above_below(j,data,aperture_array,apertures_profile_middle,aperture_peak)

    spec = []
    spec_error = []
    spec_mask = []

    for x in pix:
        ymid = aperture_array[j].trace_func(x)+offset
        profile_sigma = aperture_array[j].profile_sigma(x)
        wing = np.min([extract1d_aperture_width,3.*profile_sigma,(above0-below0)/2./2.])
#        wing = 3.
        y1 = np.long(ymid-wing)
        y2 = np.long(ymid+wing)
        sss = data[y1:y2+1,x]
        if ((wing>0.)&(len(np.where(sss.mask==False)[0])>=1)):
#            sum1 = CCDData([0.],unit=data.unit,mask=[False])
            sum1 = 0.
            sum2 = 0.
            for k in range(0,len(sss.data)):
                if sss.mask[k]==False:
                    int1 = 0.5*(scipy.special.erf(ymid/profile_sigma/np.sqrt(2.))-scipy.special.erf((ymid-(y1-0.5+k))/profile_sigma/np.sqrt(2.)))
                    int2 = 0.5*(scipy.special.erf(ymid/profile_sigma/np.sqrt(2.))-scipy.special.erf((ymid-(y1-0.5+k+1))/profile_sigma/np.sqrt(2.)))
                    weight = int2-int1
#                    sum1=sum1.add((sss[k].multiply(weight)).divide(sss.uncertainty.quantity.value[k]**2))
                    sum1 += sss.data[k]*weight/sss.uncertainty.quantity.value[k]**2
                    sum2 += weight**2/sss.uncertainty.quantity.value[k]**2
#            val = sum1.divide(sum2)
            val = sum1/sum2
            err = 1./np.sqrt(sum2)
            spec.append(val)
            spec_error.append(err)
            if ((val==val)&(err>0.)&(err==err)):
                spec_mask.append(False)
            else:
                spec_mask.append(True)
#            spec.append(sss.data[0])
#            spec_error.append(sss.uncertainty.quantity.value[0])
#            spec_mask.append(False)            
        else:
            spec.append(0.)
            spec_mask.append(True)
            spec_error.append(999.)
    spec = np.array(spec)
    spec_mask = np.array(spec_mask)
    spec_error = np.array(spec_error)

    spec_mask[np.where(((pix<aperture_array[j].trace_pixel_min)|(pix>aperture_array[j].trace_pixel_max)))[0]] = True
    return extract1d(aperture=aperture_array[j].trace_aperture,spec1d_pixel=pix,spec1d_flux=spec*data.unit,
                     spec1d_uncertainty=StdDevUncertainty(spec_error),spec1d_mask=spec_mask)


def get_apflat(data,aperture_array,apertures_profile_middle,aperture_peak,image_boundary,apmask,field_name):
    import numpy as np
    import dill as pickle
    import scipy
    from astropy.nddata import CCDData
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from astropy.nddata import StdDevUncertainty
    import matplotlib.pyplot as plt

    pix = np.arange(len(data.data[0]))
    arr2d = np.zeros((len(data.data),len(data.data[0])))
    apmask0 = ~apmask  # np.array(data.mask,dtype=bool)

    apflat_mask = mask_boundary(apmask0,image_boundary)

    for j in range(0,len(aperture_array)):
        print(j,len(aperture_array))
        if apertures_profile_middle.realvirtual[j]:
            
            above0,below0 = get_above_below(j,data,aperture_array,apertures_profile_middle,aperture_peak)

            for x in pix:
                ymid = aperture_array[j].trace_func(x)
#                print(j,x,ymid)
                profile_sigma = aperture_array[j].profile_sigma(x)
                profile_amplitude = aperture_array[j].profile_amplitude(x)
                wing = np.min([3.*profile_sigma,(above0-below0)/2./2.])
                if wing>0.:
                    y1 = np.long(ymid-wing)
                    y2 = np.long(ymid+wing)
                    if ((x<aperture_array[j].trace_pixel_min)|(x>aperture_array[j].trace_pixel_max)):
                        apflat_mask[y1:y2,x] = True
                    if ((ymid<0.)|(ymid>len(data.data))):
                        apflat_mask[y1:y2,x] = True
                    if (y1>=0)&(y2<=len(data.data)):
                        sss = data[y1:y2+1,x]
#                a.mask[y1:y2,x] = True
                        for k in range(0,len(sss.data)):
                            int1 = 0.5*(scipy.special.erf(ymid/profile_sigma/np.sqrt(2.))-scipy.special.erf((ymid-(y1-0.5+k))/profile_sigma/np.sqrt(2.)))
                            int2 = 0.5*(scipy.special.erf(ymid/profile_sigma/np.sqrt(2.))-scipy.special.erf((ymid-(y1-0.5+k+1))/profile_sigma/np.sqrt(2.)))
                            weight = int2-int1
#                            print(k,y1,y2,ymid,len(sss.data),x)
                            arr2d[y1+k,x] = weight*profile_amplitude*np.sqrt(2.*np.pi*profile_sigma**2)
#                        print(weight,profile_amplitude,profile_sigma,data.data[y1+k,x],arr2d[y1+k,x],data.data[y1+k,x]/arr2d[y1+k,x])
    residual = CCDData((data.data-arr2d)/data.uncertainty.quantity.value,
                       uncertainty=StdDevUncertainty(data.uncertainty.quantity.value/data.uncertainty.quantity.value),
                       unit=u.electron/u.electron,mask=apflat_mask)
    rms = np.sqrt(np.mean((data.data[apmask==True]-arr2d[apmask==True])**2))
    new_masked = np.where(np.abs(residual.data)>3.)  # this flags the 'ghosts'
    if len(new_masked[0])>0:
        if 'twilight' not in field_name:  # when twilight is its own flat-field spectrum, this doesn't work
            apflat_mask[new_masked] = True

#    expected_2darray = CCDData(arr2d,unit=data.unit,mask=data.mask,uncertainty=StdDevUncertainty(np.full((len(data.data),len(data.data[0])),rms,dtype='float')))
#    pickle.dump(expected_2darray,open('/nfs/nas-0-9/mgwalker.proj/m2fs/crap.pickle','wb'))
#    np.pause()
#    crap = data.uncertainty.quantity.value/arr2d
#    shite = ~apmask
    apflat = CCDData(data.data/arr2d,unit=u.electron/u.electron,mask=apflat_mask,
                     uncertainty=StdDevUncertainty(data.uncertainty.quantity.value/arr2d))

#    apflat=data.multiply(expected_2darray)    
 #   print(expected_2darray.uncertainty)
 #   print(' ')
#    apflatcorr=data.divide(apflat)
    return apflat,residual


def mask_boundary(mask0,image_boundary):
    import numpy as np

    mask = mask0
    for i in range(0,len(mask)):
        mask_low = np.long(np.interp(x=np.float(i),xp=[image_boundary.lower[q][1] for q in range(0,len(image_boundary.lower))],fp=[image_boundary.lower[q][0] for q in range(0,len(image_boundary.lower))]))
        mask_high = np.long(np.interp(x=np.float(i),xp=[image_boundary.upper[q][1] for q in range(0,len(image_boundary.upper))],fp=[image_boundary.upper[q][0] for q in range(0,len(image_boundary.upper))]))
        mask[i][0:mask_low] = True
        mask[i][mask_high:len(mask[i])] = True
    return mask


def get_above_below(j,data,aperture_array,apertures_profile_middle,aperture_peak):
    import numpy as np

    below = np.where([aperture_array[q].trace_aperture<aperture_array[j].trace_aperture for q in np.where(apertures_profile_middle.realvirtual)[0]])[0]
    above = np.where([aperture_array[q].trace_aperture>aperture_array[j].trace_aperture for q in np.where(apertures_profile_middle.realvirtual)[0]])[0]
    if len(below)>0:
        below0 = aperture_peak[np.max(below)]
    else:
        below0 = 0.
    if len(above)>0:
        above0 = aperture_peak[np.min(above)]
    else:
        above0 = np.float(len(data.data))
    return above0,below0


def get_apmask(data,aperture_array,apertures_profile_middle,aperture_peak,image_boundary):
    import numpy as np
    from astropy.nddata import CCDData
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from astropy.nddata import StdDevUncertainty

    apmask0 = np.array(data.mask,dtype=bool)
    pix = np.arange(len(apmask0[0]))

    apmask = mask_boundary(apmask0,image_boundary)

    for j in range(0,len(aperture_array)):
        if apertures_profile_middle.realvirtual[j]:

            above0,below0 = get_above_below(j,data,aperture_array,apertures_profile_middle,aperture_peak)

            for x in pix:
#                if ((x>=aperture_array[j].trace_pixel_min)&(x<=aperture_array[j].trace_pixel_max)):
                ymid = aperture_array[j].trace_func(x)
                profile_sigma = aperture_array[j].profile_sigma(x)
                wing = np.min([3.*profile_sigma,(above0-below0)/2./2.])
                if wing>0.:
                    y1 = np.long(ymid-wing)
                    y2 = np.long(ymid+wing)
                    apmask[y1:y2,x] = True
                    if ((x<aperture_array[j].trace_pixel_min)|(x>aperture_array[j].trace_pixel_max)):
                        apmask[y1:y2,x] = True
                    if ((ymid<0.)|(ymid>len(data.data))):
                        apmask[y1:y2,x] = True
    return apmask

def get_scatteredlightfunc(data,apmask,scatteredlightcorr_order,scatteredlightcorr_rejection_iterations,
                           scatteredlightcorr_rejection_sigma):
    import numpy as np
    import astropy.units as u
    from astropy.modeling import models
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    import matplotlib.pyplot as plt

    func_init = models.Polynomial2D(degree=scatteredlightcorr_order)
    fitter = fitting.LinearLSQFitter()
    y,x = np.mgrid[:len(data.data),:len(data.data[0])]

    for q in range(0,scatteredlightcorr_rejection_iterations):
        use = np.where(apmask==False)[0]
        func = fitter(func_init,x[apmask==False],y[apmask==False],data.data[apmask==False])
        rms = np.sqrt(np.mean((data.data[apmask==False]-func(x[apmask==False],y[apmask==False]))**2))
        outlier = np.where(np.abs(data.data-func(x,y))>scatteredlightcorr_rejection_sigma*rms)
#        outlier = np.where(data.data-func(x,y)>scatteredlightcorr_rejection_sigma*rms)
        apmask[outlier] = True
    return scatteredlightfunc(func=func,rms=rms)

def get_linelist(linelist_file):
    import numpy as np
    with open(linelist_file) as f:
        data = f.readlines()
    linelist_wavelength = []
    linelist_species = []
    for line in data:
        p = line.split()
        linelist_wavelength.append(float(p[0]))
        linelist_species.append(p[1])
    linelist_wavelength = np.array(linelist_wavelength)
    linelist_species = np.array(linelist_species)
    return linelist(wavelength=linelist_wavelength,species=linelist_species)


def get_wav(i,thar,extract1d_array,thar_mjd,mjd,id_lines_minlines):
    # wav array for bad apertures (all masked pixels, number of id'd lines < id_lines_minlines, etc.) will have one element [0]
    import numpy as np
    import scipy

    wav = []
    wav_rms = []
    line_sigma = []

    thar_func = []
    thar_use = []
    thar_npoints = []
    thar_rms = []
    thar_dwav = []
    thar_dpix = []
    thar_resolution = []

    for j in range(0,len(thar)):
        this = np.where(np.array([thar[j][q].aperture for q in range(0,len(thar[j]))])==extract1d_array[i].aperture)[0]
        if len(this)>0:
            if thar[j][this[0]].npoints>=id_lines_minlines:              
                thar_func.append(thar[j][this[0]].func)
                thar_npoints.append(thar[j][this[0]].npoints)
                thar_rms.append(thar[j][this[0]].rms)
                keep = np.where(thar[j][this[0]].wav.mask==False)[0]
                thar_dwav.append(np.max(thar[j][this[0]].wav)-np.min(thar[j][this[0]].wav))
                thar_dpix.append(np.max([thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep])-np.min([thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep]))
                thar_resolution.append(thar[j][this[0]].resolution)
                thar_use.append(1)
            else:
                thar_func.append(0.)
                thar_npoints.append(-999.)
                thar_rms.append(-999.)
                thar_dwav.append(-999.)
                thar_dpix.append(-999.)
                thar_resolution.append(thar[j][this[0]].resolution)
                thar_use.append(0)
    thar_use = np.array(thar_use)
    thar_npoints = np.array(thar_npoints)
    thar_rms = np.array(thar_rms)
    thar_dwav = np.array(thar_dwav)
    thar_dpix = np.array(thar_dpix)
    thar_resolution = np.array(thar_resolution)

    keep = np.where(thar_use==1)[0]
    if len(keep)>0:
        pix = extract1d_array[i].spec1d_pixel
        wav = thar_func[keep[0]](pix)#default is to use first thar function
        wav_rms = np.full(len(wav),-999.,dtype=float)
        line_sigma = thar_resolution[keep[0]](pix)*thar_dwav[keep[0]]/thar_dpix[keep[0]]
        if len(keep)>1:#if there is more than 1 thar function, then interpolate according to MJD of observation
            wav = []
            wav_rms = []
            line_sigma = []
            for j in range(0,len(pix)):
                wav0 = []
                line_sigma0 = []
                for k in range(0,len(keep)):
                    wav0.append(thar_func[keep[k]](pix[j]))
                    line_sigma0.append(thar_resolution[keep[k]](pix[j])*thar_dwav[keep[k]]/thar_dpix[keep[k]])
                wav0 = np.array(wav0)
                line_sigma0 = np.array(line_sigma0)

                def my_func(q):
                    import numpy as np
#                    print(len(thar_mjd),len(thar_rms),len(thar_use),thar_use)
                    chi2 = (wav0-(q[1]*thar_mjd[thar_use==1]+q[0]))**2/thar_rms[thar_use==1]**2
                    return np.sum(chi2)

#                        f = interpolate.interp1d(thar_mjd[thar_use],wav0,kind='linear',fill_value=(wav0[0],wav0[-1]),bounds_error=False)
                params = np.array([wav0[0],0.])
                lsq = scipy.optimize.minimize(my_func,params,method='powell',options={'maxiter':100}) # use LLSq to fit linear function to wavelength (at this pixel) as function of MJD of arc
#                        print(lsq.x,wav0,thar_mjd[thar_use],mjd,lsq.x[0]+lsq.x[1]*mjd,f(mjd))
                wav.append(lsq.x[0]+lsq.x[1]*mjd)#apply that function to get wavelength (at this pixel) at MJD of flat
                wav_rms.append(np.sqrt(np.mean((wav0-(lsq.x[0]+lsq.x[1]*mjd))**2)))
                line_sigma.append(np.mean(line_sigma0))
#                        wav.append(f(mjd))
            wav = np.array(wav)
            wav_rms = np.array(wav_rms)
            line_sigma = np.array(line_sigma)

    # line_sigma is function estimated based on stddev of gaussian fits to arc lines
    return wavcal(wav=wav,nthar=len(np.where(thar_use==1)[0]),thar_npoints=thar_npoints,thar_rms=thar_rms,wav_rms=wav_rms,mjd=mjd,line_sigma=line_sigma)


def get_throughput_continuum(twilightstack_array,twilightstack_wavcal_array,continuum_rejection_low,
                             continuum_rejection_high,continuum_rejection_iterations,continuum_rejection_order):
    import numpy as np
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative
    from copy import deepcopy
    import matplotlib.pyplot as plt
    from astropy.modeling import models
    from scipy import interpolate
    import scipy

    continuum_array = []
    wav_min = []
    wav_max = []
    npix = []
    use = []
    for i in range(0,len(twilightstack_array)):
#        print('working on aperture ',i+1,' of ',len(twilightstack_array))
        wav = twilightstack_wavcal_array[i].wav  # get_wav(i,thar,twilightstack_array,thar_mjd,mjd,id_lines_minlines)
#        print(i,len(np.where(twilightstack_array[i].spec1d_mask==False)[0]),len(wav))
        if ((len(np.where(twilightstack_array[i].spec1d_mask==False)[0])>1)&(len(wav)>0)):
            wav_min.append(np.min(wav[twilightstack_array[i].spec1d_mask==False]))
            wav_max.append(np.max(wav[twilightstack_array[i].spec1d_mask==False]))
            npix.append(len(np.where(twilightstack_array[i].spec1d_mask==False)[0]))
            use.append(1)
#            print(np.max(twilightstack_array[i].spec1d_flux[twilightstack_array[i].spec1d_mask]))
            spec1d0 = Spectrum1D(spectral_axis=wav*u.AA,flux=deepcopy(twilightstack_array[i].spec1d_flux)*u.electron,
                                 uncertainty=deepcopy(twilightstack_array[i].spec1d_uncertainty),mask=deepcopy(twilightstack_array[i].spec1d_mask))
#            plt.plot(wav[twilightstack_array[i].spec1d_mask==False],twilightstack_array[i].spec1d_flux[twilightstack_array[i].spec1d_mask==False])
#            plt.show()
#            plt.close()
            continuum0,rms0=get_continuum(spec1d0,continuum_rejection_low,continuum_rejection_high,
                                          continuum_rejection_iterations,continuum_rejection_order)
            continuum_array.append(continuum0)
        else:
            wav_min.append(-999.)
            wav_max.append(-999.)
            npix.append(0)
            use.append(0)
            continuum_array.append(models.Chebyshev1D(degree=1))

    use = np.array(use)
    wav_min = np.array(wav_min)
    wav_max = np.array(wav_max)
    npix = np.array(npix)
    continuum_array = np.array(continuum_array)
    if((len(twilightstack_array)>0)&(len(np.where(use==1)[0])>0)):
        wav0 = np.linspace(np.min(wav_min[use==1]),np.max(wav_max[use==1]),np.long(np.median(npix[use==1])))
        flux0 = []  # np.zeros(len(wav0))
        mask0 = []
        for j in range(0,len(wav0)):
            keep = np.where((wav_min<=wav0[j])&(wav_max>=wav0[j])&(use==1))[0]
            if len(keep)>0:
                flux0.append(np.median([continuum_array[q](wav0[j]) for q in keep]))
                mask0.append(False)
            else:
                flux0.append(0.)
                mask0.append(True)
        flux0 = np.array(flux0)    
        mask0 = np.array(mask0)
        spec1d0 = Spectrum1D(spectral_axis=wav0*u.AA,flux=flux0*u.electron,mask=mask0)
        continuum0,rms0 = get_continuum(spec1d0,continuum_rejection_low,continuum_rejection_high,
                                        continuum_rejection_iterations,continuum_rejection_order)
        
#        plt.plot(spec1d0.spectral_axis,spec1d0.flux)
#        plt.plot(spec1d0.spectral_axis,continuum0(spec1d0.spectral_axis.value))
#        plt.show()
#        plt.close()
#        np.pause()
    else:
        continuum0 = models.Chebyshev1D(degree=1)
        rms0 = -999.

    return continuum_array,continuum0


def get_throughputcorr(i,throughput_continuum,extract1d_array,twilightstack_array,
                       twilightstack_continuum_array,wavcal_array):
    # thar,id_lines_minlines,mjd,thar_mjd,flat_mjd):
    import numpy as np
    from specutils.spectra import Spectrum1D
    import astropy.units as u
    from specutils.fitting import find_lines_threshold
    from specutils.fitting import find_lines_derivative
    from copy import deepcopy
    import matplotlib.pyplot as plt
    from astropy.modeling import models
    from scipy import interpolate
    import scipy

    continuum_array = []
    wav_min = []
    wav_max = []
    npix = []

    spec1d0 = Spectrum1D(spectral_axis=extract1d_array[i].spec1d_pixel*u.AA,flux=extract1d_array[i].spec1d_flux,
                         uncertainty=extract1d_array[i].spec1d_uncertainty,mask=extract1d_array[i].spec1d_mask)
    throughput = np.full((len(spec1d0.flux.value)),1.,dtype='float')*u.electron/u.electron
    throughputcorr = deepcopy(extract1d_array[i])
    throughputcorr.spec1d_mask = np.full(len(throughputcorr.spec1d_flux),True,dtype=bool)

    if len(twilightstack_array) != len(wavcal_array):
        print(' twilightstack_array different length than twilightstack_wavcal_array!!!')
        np.pause()
#    for j in range(0,len(twilightstack_array)):
#        twilightstack_array[j].spec1d_flux=twilightstack_continuum_array[j](twilightstack_wavcal_array[j].wav)*u.electron

    if len(np.where(np.array(extract1d_array[i].spec1d_mask)==False)[0])>0:
        this1 = np.where(np.array([twilightstack_array[q].aperture for q in range(0,len(twilightstack_array))])==extract1d_array[i].aperture)[0]
        if len(this1)>0:
            throughput = twilightstack_continuum_array[this1[0]](wavcal_array[this1[0]].wav)/throughput_continuum(wavcal_array[this1[0]].wav)*u.electron/u.electron
#                throughput00 = np.interp(wav0,twilightstack_wavcal_array[this1[0]].wav,throughput0)
#                throughput = np.interp(wavcal_array[i].wav,wav0,throughput00)*u.electron/u.electron
            throughputcorr0 = spec1d0.divide(throughput)
            throughputcorr = extract1d(aperture=extract1d_array[i].aperture,spec1d_pixel=extract1d_array[i].spec1d_pixel,
                                       spec1d_flux=throughputcorr0.flux,spec1d_uncertainty=throughputcorr0.uncertainty,spec1d_mask=throughputcorr0.mask)

#    else:
#        print('')

    return throughput,throughputcorr


def get_plugmap(header,throughputcorr_array,fibermap,ccd,fiber_changes):
    # for older hecto frames, interprets iraf headers to get wavelength solutions for each aperture
    from astropy.io import fits
    import numpy as np
    import re
    from astropy import time, coordinates as coord, units as u
    from astropy.coordinates import SkyCoord, EarthLocation

    changed = np.zeros(len(fibermap),dtype=bool)
    if not('none' in fiber_changes):
        for i in range(0,len(fiber_changes)):
            old,new = fiber_changes[i].split('>')
            print('implementing fiber plugging change ',fiber_changes[i])
            print('notation: [fiber being replaced > fiber doing the replacing]')
            new_j = -1
            old_j = -1
            for j in range(0,len(fibermap)):
                if ((old in fibermap[j])&(not('dead' in fibermap[j]))):
                    old_j = j
                if ((new in fibermap[j])&(not('dead' in fibermap[j]))):
                    new_j = j
            if old_j<0:
                print('ERROR--DID NOT FIND REPLACED FIBER IN FIBERMAP')
                np.pause()
            if new_j>=0:
                if new!='unplugged':
                    fibermap[new_j] = new+fibermap[old_j].split(old)[1]
            fibermap[old_j] = old+' unplugged - - - I - - - - - - - -'

    if ccd=='r':
        channel0 = 'R'
    if ccd=='b':
        channel0 = 'B'
    resolution0 = header['slide']
    filter0 = header['filter']

    aperture = []
    radeg = []
    decdeg = []
    objtype = []
    expid = []
    icode = []
    rcode = []
    xfocal = []
    yfocal = []
    channel = []
    resolution = []
    filt = []
    channelcassfib = []

    for i in range(0,len(throughputcorr_array)):
        ap = 129-throughputcorr_array[i].aperture
        if channel0=='B':
            if ((ap>=1)&(ap<=16)):
                cass = '8'
                fib = str(ap).zfill(2)
            if ((ap>=17)&(ap<=32)):
                cass = '7'
                fib = str(ap-16).zfill(2)
            if ((ap>=33)&(ap<=48)):
                cass = '6'
                fib = str(ap-32).zfill(2)
            if ((ap>=49)&(ap<=64)):
                cass = '5'
                fib = str(ap-48).zfill(2)
            if ((ap>=65)&(ap<=80)):
                cass = '4'
                fib = str(ap-64).zfill(2)
            if ((ap>=81)&(ap<=96)):
                cass = '3'
                fib = str(ap-80).zfill(2)
            if ((ap>=97)&(ap<=112)):
                cass = '2'
                fib = str(ap-96).zfill(2)
            if ((ap>=113)&(ap<=128)):
                cass = '1'
                fib = str(ap-112).zfill(2)
        if channel0=='R':
            if ((ap>=1)&(ap<=16)):
                cass = '1'
                fib = str(ap).zfill(2)
            if ((ap>=17)&(ap<=32)):
                cass = '2'
                fib = str(ap-16).zfill(2)
            if ((ap>=33)&(ap<=48)):
                cass = '3'
                fib = str(ap-32).zfill(2)
            if ((ap>=49)&(ap<=64)):
                cass = '4'
                fib = str(ap-48).zfill(2)
            if ((ap>=65)&(ap<=80)):
                cass = '5'
                fib = str(ap-64).zfill(2)
            if ((ap>=81)&(ap<=96)):
                cass = '6'
                fib = str(ap-80).zfill(2)
            if ((ap>=97)&(ap<=112)):
                cass = '7'
                fib = str(ap-96).zfill(2)
            if ((ap>=113)&(ap<=128)):
                cass = '8'
                fib = str(ap-112).zfill(2)

        channel_cass_fib = channel0+cass+'-'+fib
        fibermap_item0 = [q for q in fibermap if ((channel_cass_fib in q)&(not('deadfibers' in q)))][0]
        fibermap_item = fibermap_item0.split()
        coords_string = fibermap_item[2]+' '+fibermap_item[3]
#        coords_string = header[cassfib].replace('_',' ')
        
        t = -999  # for consistency with hecto fits format
        r = -999  # for consistency with hecto fits format
        x = -999  # for consistency with hecto fits format
        y = -999  # for consistency with hecto fits format
        obj = coords_string
        if ((fibermap_item[5]=='T')|(fibermap_item[5]=='O')):
#        if ':' in coords_string:
            c = SkyCoord(coords_string,unit=(u.hourangle,u.deg))
            aperture.append(np.long(throughputcorr_array[i].aperture))
            radeg.append(c.ra.degree)
            decdeg.append(c.dec.degree)
            expid.append(obj)
            rcode.append(r)
            xfocal.append(x)
            yfocal.append(y)
            objtype0 = 'TARGET'
#            if 'sky' in obj: objtype0='SKY'
#            if 'unused' in obj: objtype0='UNUSED'
#            if 'reject' in obj: objtype0='REJECTED'
#            if 'SKY' in obj: objtype0='SKY'
#            if 'UNUSED' in obj: objtype0='UNUSED'
#            if 'REJECT' in obj: objtype0='REJECTED'
#            if objtype0=='SKY': 
#                icode.append(-1)
#            if objtype0=='TARGET': 
            icode.append(1)
#            if objtype0=='UNUSED':
#                icode.append(0)
            objtype.append(objtype0)
            channel.append(channel0)
            resolution.append(resolution0)
            filt.append(filter0)
            channelcassfib.append(channel_cass_fib)

        elif fibermap_item[5]=='S':
#        if (('sky' in coords_string)|('Sky' in coords_string)):
            c = SkyCoord(coords_string,unit=(u.hourangle,u.deg))
            aperture.append(np.long(throughputcorr_array[i].aperture))
            radeg.append(c.ra.degree)
            decdeg.append(c.dec.degree)
            expid.append(obj)
            rcode.append(r)
            xfocal.append(x)
            yfocal.append(y)
            objtype0 = 'SKY'
#            if 'sky' in obj: objtype0='SKY'
#            if 'unused' in obj: objtype0='UNUSED'
#            if 'reject' in obj: objtype0='REJECTED'
#            if 'SKY' in obj: objtype0='SKY'
#            if 'UNUSED' in obj: objtype0='UNUSED'
#            if 'REJECT' in obj: objtype0='REJECTED'
#            if objtype0=='SKY': 
#            icode.append(-1)
#            if objtype0=='TARGET': 
#            icode.append(1)
#            if objtype0=='UNUSED':
            icode.append(0)
            objtype.append(objtype0)
            channel.append(channel0)
            resolution.append(resolution0)
            filt.append(filter0)
            channelcassfib.append(channel_cass_fib)
        else:
#            c = SkyCoord(coords_string,unit=(u.hourangle,u.deg))
            aperture.append(np.long(throughputcorr_array[i].aperture))
            radeg.append(-999)
            decdeg.append(-999)
            expid.append(obj)
            rcode.append(r)
            xfocal.append(x)
            yfocal.append(y)
            objtype0 = 'unused'
#            if 'sky' in obj: objtype0='SKY'
#            if 'unused' in obj: objtype0='UNUSED'
#            if 'reject' in obj: objtype0='REJECTED'
#            if 'SKY' in obj: objtype0='SKY'
#            if 'UNUSED' in obj: objtype0='UNUSED'
#            if 'REJECT' in obj: objtype0='REJECTED'
#            if objtype0=='SKY': 
            icode.append(-1)
#            if objtype0=='TARGET': 
#            icode.append(1)
#            if objtype0=='UNUSED':
#                icode.append(0)
            objtype.append(objtype0)
            channel.append(channel0)
            resolution.append(resolution0)
            filt.append(filter0)
            channelcassfib.append(channel_cass_fib)

    icode = np.array(icode,dtype='int')
    objtype = np.array(objtype)
    expid = np.array(expid)
    radeg = np.array(radeg)
    decdeg = np.array(decdeg)
    aperture = np.array(aperture)
    rcode = np.array(rcode)
    xfocal = np.array(xfocal)
    yfocal = np.array(yfocal)
    channel = np.array(channel)
    resolution = np.array(resolution)
    filt = np.array(filt)

    bcode = objtype  # redundancy seems to be built in at CfA for some reason
    rmag = np.zeros(len(aperture))
    rapmag = np.zeros(len(aperture))
    frames = np.zeros(len(aperture),dtype='int')
    mag = np.zeros((len(aperture),5))

    cols = fits.ColDefs([fits.Column(name='EXPID',format='A100',array=expid),
                         fits.Column(name='OBJTYPE',format='A6',array=objtype),
                         fits.Column(name='RA',format='D',array=radeg),
                         fits.Column(name='DEC',format='D',array=decdeg),
                         fits.Column(name='APERTURE',format='I',array=aperture),
                         fits.Column(name='RMAG',format='D',array=rmag),
                         fits.Column(name='RAPMAG',format='D',array=rapmag),
                         fits.Column(name='ICODE',format='D',array=icode),
                         fits.Column(name='RCODE',format='D',array=rcode),
                         fits.Column(name='BCODE',format='A6',array=bcode),
                         fits.Column(name='MAG',format='5D',array=mag),
                         fits.Column(name='XFOCAL',format='D',array=xfocal),
                         fits.Column(name='YFOCAL',format='D',array=yfocal),
                         fits.Column(name='FRAMES',format='B',array=frames),
                         fits.Column(name='CHANNEL',format='A100',array=channel),
                         fits.Column(name='RESOLUTION',format='A100',array=resolution),
                         fits.Column(name='FILTER',format='A100',array=filt),
                         fits.Column(name='CHANNEL_CASSETTE_FIBER',format='A100',array=channelcassfib)])

    plugmap_table_hdu = fits.FITS_rec.from_columns(cols)
    return plugmap_table_hdu


def get_meansky(throughputcorr_array,wavcal_array,plugmap):
    import numpy as np
    import scipy
    import astropy.units as u
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    import matplotlib.pyplot as plt
    from astropy.nddata import StdDevUncertainty

    wav_min = []
    wav_max = []
    npix = []
    use = []
    for i in range(0,len(throughputcorr_array)):
        pix = throughputcorr_array[i].spec1d_pixel
        wav = np.array(wavcal_array[i].wav)
        keep = np.where(throughputcorr_array[i].spec1d_mask==False)[0]
        if ((len(wav)>0)&(len(keep)>0)):
            wav_min.append(np.min(wav[throughputcorr_array[i].spec1d_mask==False]))
            wav_max.append(np.max(wav[throughputcorr_array[i].spec1d_mask==False]))
            npix.append(len(np.where(throughputcorr_array[i].spec1d_mask==False)[0]))
            use.append(1)
        else:
            wav_min.append(1.e+10)
            wav_max.append(-1.e+10)
            npix.append(0)
            use.append(0)
    wav_min = np.array(wav_min)
    wav_max = np.array(wav_max)
    npix = np.array(npix)
    use = np.array(use)

    meansky = Spectrum1D(spectral_axis=np.array([0.])*u.AA,flux=np.array([0.])*u.electron,
                         uncertainty=StdDevUncertainty(np.array([np.inf])),mask=np.array([True]))
    if len(throughputcorr_array)>0:
        if len(wav_min[use==1])>0:
            wav0 = np.linspace(np.min(wav_min[use==1]),np.max(wav_max[use==1]),np.long(np.median(npix[use==1])*10))
#            id_nlines = np.array([len(np.where(id_lines_array[q].wav>0.)[0]) for q in range(0,len(id_lines_array))],dtype='int')
            skies = np.where(((plugmap['objtype']=='SKY')|(plugmap['objtype']=='unused'))&(np.array([len(wavcal_array[q].wav) for q in range(0,len(wavcal_array))])>0))[0]
            targets = np.where((plugmap['objtype']=='TARGET')&(np.array([len(wavcal_array[q].wav) for q in range(0,len(wavcal_array))])>0))[0]
            
            sky_flux_array = []
            sky_err_array = []
            sky_mask_array = []

            for i in range(0,len(skies)):
#            use = np.where(throughputcorr_array[skies[i]].spec1d_mask==False)[0]
                wav = wavcal_array[skies[i]].wav
                if len(wav)>0:
                    sky_flux_array.append(np.interp(wav0,wav,throughputcorr_array[skies[i]].spec1d_flux))
                    sky_err_array.append(np.interp(wav0,wav,throughputcorr_array[skies[i]].spec1d_uncertainty.quantity.value))
                    sky_mask_array.append(np.interp(wav0,wav,throughputcorr_array[skies[i]].spec1d_mask))
            sky_flux_array = np.array(sky_flux_array)
            sky_err_array = np.array(sky_err_array)
            sky_mask_array = np.array(sky_mask_array)

            sky0_flux = []
            sky0_err = []
            sky0_mask = []
            for i in range(0,len(wav0)):
                vec = np.array([sky_flux_array[q][i] for q in range(0,len(sky_flux_array))],dtype='float')
                vec_mask = np.array([sky_mask_array[q][i] for q in range(0,len(sky_flux_array))],dtype='bool')
                for j in range(0,len(vec)):
                    if wav_min[skies[j]]>wav0[i]:
                        vec_mask[j] = True
                    if wav_max[skies[j]]<=wav0[i]:
                        vec_mask[j] = True
                med = np.median(vec[vec_mask==False])
                mad = np.median(np.abs(vec[vec_mask==False]-med))*1.4826*np.sqrt(np.pi/2.)/np.sqrt(np.float(len(np.where(vec_mask==False)[0])))

                sky0_flux.append(med)
                sky0_err.append(mad)
                if med!=med:
                    sky0_mask.append(True)
                else:
                    sky0_mask.append(False)
            sky0_flux = np.array(sky0_flux)
            sky0_err = np.array(sky0_err)
            sky0_mask = np.array(sky0_mask)
        else:
            wav0 = [0.]
            sky0_flux = np.array([0.])
            sky0_err = np.array([0.])
            sky0_mask = np.array([True])

        meansky = Spectrum1D(spectral_axis=wav0*u.AA,flux=sky0_flux*u.electron,
                             uncertainty=StdDevUncertainty(sky0_err),mask=sky0_mask)
    return meansky


def get_skysubtract(meansky,i,throughputcorr_array,wavcal_array):
    import numpy as np
    import scipy
    import astropy.units as u
    from specutils.spectra import Spectrum1D
    from astropy.modeling import models,fitting
    import matplotlib.pyplot as plt
    from astropy.nddata import StdDevUncertainty
    from copy import deepcopy

    spec1d0 = Spectrum1D(spectral_axis=throughputcorr_array[i].spec1d_pixel*u.AA,
                         flux=throughputcorr_array[i].spec1d_flux,
                         uncertainty=throughputcorr_array[i].spec1d_uncertainty,
                         mask=throughputcorr_array[i].spec1d_mask)
    sky0 = Spectrum1D(spectral_axis=throughputcorr_array[i].spec1d_pixel*u.AA,
                      flux=np.zeros(len(throughputcorr_array[i].spec1d_mask))*u.electron,
                      uncertainty=np.full(len(throughputcorr_array[i].spec1d_mask),np.inf),
                      mask=np.full(len(throughputcorr_array[i].spec1d_mask),True))
    wav0 = meansky.spectral_axis.value
    wav = wavcal_array[i].wav
    skysubtract0 = deepcopy(spec1d0)
    skysubtract0.mask = np.full(len(skysubtract0.flux),True,dtype=bool)

    if len(wav)>0:
        sky_flux = np.interp(wav,wav0,meansky.flux)
        sky_err = np.interp(wav,wav0,meansky.uncertainty.quantity.value)
        sky_mask = np.interp(wav,wav0,meansky.mask)
        sky0 = Spectrum1D(spectral_axis=throughputcorr_array[i].spec1d_pixel*u.AA,
                          flux=sky_flux,uncertainty=StdDevUncertainty(sky_err),mask=sky_mask)
        skysubtract0 = spec1d0.subtract(sky0)
    sky = extract1d(throughputcorr_array[i].aperture,spec1d_pixel=throughputcorr_array[i].spec1d_pixel,
                    spec1d_flux=sky0.flux,spec1d_uncertainty=sky0.uncertainty,spec1d_mask=sky0.mask)
    skysubtract = extract1d(throughputcorr_array[i].aperture,spec1d_pixel=throughputcorr_array[i].spec1d_pixel,
                            spec1d_flux=skysubtract0.flux,spec1d_uncertainty=skysubtract0.uncertainty,
                            spec1d_mask=skysubtract0.mask)
    return sky,skysubtract


def get_thars(q,thar,temperature):
    import numpy as np

    wav = []
    pix = []
    wav_func = []
    vel_func = []
    pix_std = []
    wav_func_std = []
    vel_func_std = []

    ap = thar[0][q].aperture
    keep_ap = True
    for xxx in range(0,len(thar)):
        if not thar[0][q].aperture in [thar[xxx][jjj].aperture for jjj in range(0,len(thar[xxx]))]:
            keep_ap = False
        this_ap = np.where(np.array([thar[xxx][jjj].aperture for jjj in range(0,len(thar[xxx]))])==thar[0][q].aperture)[0][0]
        if len(thar[xxx][this_ap].wav)<=1:
            keep_ap = False
    if keep_ap:

        this1 = np.where(np.array([thar[0][jjj].aperture for jjj in range(0,len(thar[0]))])==ap)[0][0]
        for yyy in range(0,len(thar[0][this1].wav)):
#            print(thar[0][this1].wav.mask)
            if len(thar[0][this1].wav)>1:
                if thar[0][this1].wav.mask[yyy]==False:
                    pix00 = thar[0][this1].fit_lines.fit[yyy].mean.value
                    wav_func00 = thar[0][this1].func(pix00)
                    keep_wav = True
                    pix0 = []
                    wav_func0 = []
                    vel_func0 = []
                    for zzz in range(0,len(thar)):
                        this2 = np.where(np.array([thar[zzz][jjj].aperture for jjj in range(0,len(thar[zzz]))])==ap)[0][0]
                        if not thar[0][this1].wav[yyy] in thar[zzz][this2].wav:
                            keep_wav = False
                        else:
                            this_wav = np.where(thar[zzz][this2].wav==thar[0][this1].wav[yyy])[0][0]
                            pix0.append(thar[zzz][this2].fit_lines.fit[this_wav].mean.value)
                            wav_func0.append(thar[zzz][this2].func(pix00))
                            if wav_func00==0.:
                                vel_func0.append(np.inf)
                            else:
                                vel_func0.append((wav_func0[len(wav_func0)-1]-wav_func00)/wav_func00*3.e+5)
                    pix0 = np.array(pix0)
                    wav_func0 = np.array(wav_func0)
                    vel_func0 = np.array(vel_func0)
                    if keep_wav:
                        wav.append(thar[0][this1].wav[yyy])
                        pix.append(pix0)
                        wav_func.append(wav_func0)
                        vel_func.append(vel_func0)
#                                pix.append(thar[][].fit_lines.fit[].mean.value)
#                                wav_func.append(thar[][]/func(pix[0]))
#                                vel_func.append((wav_func[len(wav_func)-1]-wav_func[0])/wav_func[0]*3.e+5)
        wav = np.array(wav)
        pix = np.array(pix)
        wav_func = np.array(wav_func)
        vel_func = np.array(vel_func)
        for qqq in range(0,len(wav)):
            pix_std.append(np.std(pix[qqq]))
            wav_func_std.append(np.std(wav_func[qqq]))
            vel_func_std.append(np.std(vel_func[qqq]))
        pix_std = np.array(pix_std)
        wav_func_std = np.array(wav_func_std)
        vel_func_std = np.array(vel_func_std)
        return thars(aperture=ap,wav=wav,pix=pix,wav_func=wav_func,vel_func=vel_func,pix_std=pix_std,
                     wav_func_std=wav_func_std,vel_func_std=vel_func_std,temperature=temperature)

    
def weightedmeanspec(stack00):
    import numpy as np

    spec = []
    spec_err = []
    spec_mask = []

    for i in range(0,len(stack00[0].data)):
        keep = np.where([stack00[q].mask[i]==False for q in range(0,len(stack00))])[0]
        if len(keep)>0:
            sum1 = np.sum(np.array([stack00[q].data[i] for q in keep])/np.array([stack00[q].uncertainty.quantity.value[i] for q in keep])**2)
            sum2 = np.sum(1./np.array([stack00[q].uncertainty.quantity.value[i] for q in keep])**2)
            spec.append(sum1/sum2)
            spec_err.append(np.sqrt(1./sum2))
            spec_mask.append(False)
#            sum1 = np.sum(np.array([stack00[q].data[i] for q in keep]))
#            sum2 = np.sum(np.array([stack00[q].uncertainty.quantity.value[i] for q in keep])**2)
#            spec.append(sum1)
#            spec_err.append(np.sqrt(sum2))
#            spec_mask.append(False)
        else:
            spec.append(0.)
            spec_err.append(-999.)
            spec_mask.append(True)
    spec = np.array(spec)
    spec_err = np.array(spec_err)
    spec_mask = np.array(spec_mask)
    return spec,spec_err,spec_mask


def get_stack(stack0,j):
    from astropy.nddata import CCDData
    from specutils.spectra import Spectrum1D
    from astropy import units as u
    from astropy.nddata import StdDevUncertainty
    import numpy as np

    stack00 = []
    for q in range(0,len(stack0)):
#        print(len(stack0),len(stack0[0]),len(stack0[q]),q,j)
        stack00.append(CCDData(stack0[q][j].spec1d_flux,uncertainty=stack0[q][j].spec1d_uncertainty,mask=stack0[q][j].spec1d_mask))
#        print(np.where(stack0[q][j].spec1d_mask==False)[0])
    spec,spec_err,spec_mask = weightedmeanspec(stack00)
    return extract1d(aperture=stack0[q][j].aperture,spec1d_pixel=stack0[q][j].spec1d_pixel*u.AA,spec1d_flux=spec*u.electron,
                     spec1d_uncertainty=StdDevUncertainty(spec_err),spec1d_mask=spec_mask)


def get_thar(datadir,utdate,ccd,tharfile,hires_exptime,medres_exptime,field_name,use_flat):
    import astropy
    import dill as pickle
    from astropy import time, coordinates as coord, units as u
    from astropy.coordinates import SkyCoord, EarthLocation
    import numpy as np

    lco = coord.EarthLocation.from_geodetic(lon=-70.6919444*u.degree,lat=-29.0158333*u.degree,height=2380.*u.meter)
    thar = []
    lines = []
    temperature = []
    thar_exptime = []
    thar_mjd = []
    for j in range(0,len(tharfile)):
        root0 = datadir+utdate+'/'+ccd+str(tharfile[j]).zfill(4)
        data_file = root0+'_stitched.fits'
        if use_flat:
            id_lines_array_file = root0+'_id_lines_array.pickle'
        else:
            id_lines_array_file = root0+'_id_lines_array_noflat.pickle'
        id_lines_array = pickle.load(open(id_lines_array_file,'rb'))
        data = astropy.nddata.CCDData.read(data_file)  # format is such that data.data[:,0] has column-0 value in all rows
        time0 = [data.header['DATE-OBS']+'T'+data.header['UT-TIME'],data.header['DATE-OBS']+'T'+data.header['UT-END']]
        times = time.Time(time0,location=lco,precision=6)
        filtername = data.header['FILTER']
#        print(data.header['filter'],data.header['exptime'])
        if filtername=='Mgb_Rev2':
            filtername = 'Mgb_HiRes'
        if filtername=='Mgb_HiRes':
            if ((np.float(data.header['exptime'])>=hires_exptime)|('twilight' in field_name)):
                thar_mjd.append(np.mean(times.mjd))
                thar.append(id_lines_array)
#                lines.append([id_lines_array[q].wav[id_lines_array[q].wav.mask==False] for q in range(0,len(id_lines_array))])
                temperature.append(data.header['T-DOME'])
                thar_exptime.append(np.float(data.header['exptime']))
        if filtername=='Mgb_MedRes':
            if ((np.float(data.header['exptime'])>=medres_exptime)|('twilight' in field_name)):
                thar_mjd.append(np.mean(times.mjd))
                thar.append(id_lines_array)
#                lines.append([id_lines_array[q].wav.data[id_lines_array[q].wav.mask==False] for q in range(0,len(id_lines_array))])
                temperature.append(data.header['T-DOME'])
                thar_exptime.append(np.float(data.header['exptime']))

    temperature = np.array(temperature)
    thar_exptime = np.array(thar_exptime)
    thar_mjd = np.array(thar_mjd)

    if len(thar)==0:
        print('ERROR: no qualifying ThArNe exposures for this resolution!!!!!')
        np.pause()
    if np.std(thar_exptime)>0.:
        if filtername=='Mgb_MedRes':
            keep = np.where(thar_exptime==np.min(thar_exptime))[0]
            thar0 = []
            for qqq in keep:
                thar0.append(thar[qqq])
            thar = thar0
            thar_mjd = thar_mjd[keep]
            thar_exptime = thar_exptime[keep]

    return(thar,lines,temperature,thar_exptime,thar_mjd)
    


def get_hdul(data,skysubtract_array,sky_array,wavcal_array,plugmap,m2fsrun,field_name,thar,temperature,aperture_array,linelist):
    import numpy as np

    from astropy.io import fits
    from astropy import time, coordinates as coord, units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    from copy import deepcopy

    data_array = []
    var_array = []
    mask_array = []
    sky = []
    wav_array = []
    mjd_array = []
    hjd_array = []
    vheliocorr_array = []
    snratio_array = []
    thar_npoints_array = []
    thar_rms_array = []
    thar_resolution_array = []
    thar_wav_min_array = []
    thar_wav_max_array = []
    temperature_array = []
    row_array = []

    temperature0 = []
    for i in range(0,len(temperature)):
        temperature0.append(float(temperature[i]))

    for i in range(0,len(plugmap)):
        temperature_array.append(temperature0)
        thar_npoints = []
        thar_rms = []
        thar_resolution = []
        thar_wav_min = []
        thar_wav_max = []
        for j in range(0,len(thar)):
            this = np.where([thar[j][q].aperture for q in range(0,len(thar[j]))]==plugmap['aperture'][i])[0]
            if len(this)>0:
                thar_npoints.append(thar[j][this[0]].npoints)
                thar_rms.append(thar[j][this[0]].rms)
                thar_resolution0 = np.float(-999.)
                thar_wav_min0 = np.float(-999.)
                thar_wav_max0 = np.float(-999.)
                if len(thar[j][this[0]].wav)>1:
                    keep = np.where(thar[j][this[0]].wav.mask==False)[0]
                    if len(keep)>1:
                        dwav = np.max(thar[j][this[0]].wav[keep])-np.min(thar[j][this[0]].wav[keep])
                        dpix = np.max([thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep])-np.min([thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep])
                        central_wavelength = np.min(linelist.wavelength)+(np.max(linelist.wavelength)-np.min(linelist.wavelength))/2.
#                        x = [thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep]
                        thar_resolution0 = thar[j][this[0]].resolution(np.interp(x=central_wavelength,xp=[thar[j][this[0]].wav[q] for q in keep],fp=[thar[j][this[0]].fit_lines.fit[q].mean.value for q in keep]))*dwav/dpix
                        thar_wav_min0 = np.min(np.ma.compressed(thar[j][this[0]].wav))
                        thar_wav_max0 = np.max(np.ma.compressed(thar[j][this[0]].wav))

                thar_resolution.append(thar_resolution0)
                thar_wav_min.append(thar_wav_min0)
                thar_wav_max.append(thar_wav_max0)

        thar_npoints_array.append(thar_npoints)
        thar_rms_array.append(thar_rms)
        thar_resolution_array.append(thar_resolution)
        thar_wav_min_array.append(thar_wav_min)
        thar_wav_max_array.append(thar_wav_max)

        if plugmap['objtype'][i]=='TARGET':
            coords = coord.SkyCoord(plugmap['ra'][i],plugmap['dec'][i],unit=(u.deg,u.deg),frame='icrs')
            lco = coord.EarthLocation.from_geodetic(lon=-70.6919444*u.degree,lat=-29.0158333*u.degree,height=2380.*u.meter)
            times = time.Time(wavcal_array[0].mjd,format='mjd',location=lco)
            light_travel_time_helio = times.light_travel_time(coords,'heliocentric')
            mjd_array.append(wavcal_array[0].mjd)
            hjd_array.append((Time(wavcal_array[0].mjd,format='mjd').jd+light_travel_time_helio).value)
            vheliocorr_array.append((coords.radial_velocity_correction('heliocentric',obstime=times).to(u.km/u.s)).value)
        else:
            mjd_array.append(wavcal_array[0].mjd)
            hjd_array.append(wavcal_array[0].mjd+2400000.)
            vheliocorr_array.append(0.)
        this = np.where([skysubtract_array[q].aperture for q in range(0,len(skysubtract_array))]==plugmap['aperture'][i])[0][0]
        snratio_array.append(np.median(skysubtract_array[this].spec1d_flux.value[skysubtract_array[this].spec1d_mask==False]/skysubtract_array[this].spec1d_uncertainty.quantity.value[skysubtract_array[this].spec1d_mask==False]))
        this = np.where([aperture_array[q].trace_aperture for q in range(0,len(aperture_array))]==plugmap['aperture'][i])[0][0]
        row = aperture_array[this].trace_func(np.float(len(data.data[0])/2.))
        row_array.append(row)
    snratio_array = np.array(snratio_array)
    m2fsrun_array = np.full(len(snratio_array),m2fsrun,dtype='a100')
    field_name_array = np.full(len(snratio_array),field_name,dtype='a100')
    thar_npoints_array = np.array(thar_npoints_array,dtype=np.object_)
    thar_rms_array = np.array(thar_rms_array,dtype=np.object_)
    thar_resolution_array = np.array(thar_resolution_array,dtype=np.object_)
    thar_wav_min_array = np.array(thar_wav_min_array,dtype=np.object_)
    thar_wav_max_array = np.array(thar_wav_max_array,dtype=np.object_)
    temperature_array = np.array(temperature_array,dtype=np.object_)
    row_array = np.array(row_array)


    cols = fits.ColDefs([fits.Column(name='EXPID',format='A100',array=plugmap['expid']),
                         fits.Column(name='OBJTYPE',format='A6',array=plugmap['objtype']),
                         fits.Column(name='RA',format='D',array=plugmap['ra']),
                         fits.Column(name='DEC',format='D',array=plugmap['dec']),
                         fits.Column(name='APERTURE',format='I',array=plugmap['aperture']),
                         fits.Column(name='RMAG',format='D',array=plugmap['rmag']),
                         fits.Column(name='RAPMAG',format='D',array=plugmap['rapmag']),
                         fits.Column(name='ICODE',format='D',array=plugmap['icode']),
                         fits.Column(name='RCODE',format='D',array=plugmap['rcode']),
                         fits.Column(name='BCODE',format='A6',array=plugmap['bcode']),
                         fits.Column(name='MAG',format='5D',array=plugmap['mag']),
                         fits.Column(name='XFOCAL',format='D',array=plugmap['xfocal']),
                         fits.Column(name='YFOCAL',format='D',array=plugmap['yfocal']),
                         fits.Column(name='FRAMES',format='B',array=plugmap['frames']),
                         fits.Column(name='CHANNEL',format='A100',array=plugmap['channel']),
                         fits.Column(name='RESOLUTION',format='A100',array=plugmap['resolution']),
                         fits.Column(name='FILTER',format='A100',array=plugmap['filter']),
                         fits.Column(name='CHANNEL_CASSETTE_FIBER',format='A100',array=plugmap['channel_cassette_fiber']),
                         fits.Column(name='MJD',format='D',array=mjd_array),
                         fits.Column(name='HJD',format='D',array=hjd_array),
                         fits.Column(name='vheliocorr',format='d',array=vheliocorr_array),
                         fits.Column(name='SNRATIO',format='d',array=snratio_array),
                         fits.Column(name='run_id',format='A100',array=m2fsrun_array),
                         fits.Column(name='field_name',format='A100',array=field_name_array),
                         fits.Column(name='wav_npoints',format='PI()',array=thar_npoints_array),
                         fits.Column(name='wav_rms',format='PD()',array=thar_rms_array),
                         fits.Column(name='wav_resolution',format='PD()',array=thar_resolution_array),
                         fits.Column(name='wav_min',format='PD()',array=thar_wav_min_array),
                         fits.Column(name='wav_max',format='PD()',array=thar_wav_max_array),
                         fits.Column(name='row',format='D',array=row_array),
                         fits.Column(name='temperature',format='PD()',array=temperature_array)])
    table_hdu = fits.FITS_rec.from_columns(cols)

    if len(skysubtract_array)>0:
        for i in range(0,len(plugmap)):
            this = np.where(np.array([skysubtract_array[q].aperture for q in range(0,len(skysubtract_array))])==plugmap['aperture'][i])[0]
            if len(this)>0:
                mask0 = deepcopy(skysubtract_array[this[0]].spec1d_mask)
                if len(np.where(skysubtract_array[this[0]].spec1d_mask==False)[0])>0:   # mask pixels outside valid extrapolation region of wavelength solution
                    best = np.where(thar_npoints_array[i]==np.max(thar_npoints_array[i]))[0]
                    spec_range = thar_wav_max_array[i][best]-thar_wav_min_array[i][best]
                    best2 = np.where(spec_range==np.min(spec_range))[0][0]
                    if thar_npoints_array[i][best[best2]]==0:
                        print(thar_npoints_array[i],this,wavcal_array[this[0]].wav)
                        np.pause()
#                    if thar_npoints_array[i][best[best2]]>0:###why are there unmasked cases where thar_npoints=0 and/or wavcal.wav=[]
#                        if len(wavcal_array[this[0]].wav)>0:
                    extra = (thar_wav_max_array[i][best[best2]]-thar_wav_min_array[i][best[best2]])/np.float(thar_npoints_array[i][best[best2]])
                    lambdamin = thar_wav_min_array[i][best[best2]]-extra
                    lambdamax = thar_wav_max_array[i][best[best2]]+extra
#                            print(wavcal_array[this[0]].wav)
                    new_mask = np.where((wavcal_array[this[0]].wav<lambdamin)|(wavcal_array[this[0]].wav>lambdamax))[0]
                    mask0[new_mask] = True
                data_array.append(skysubtract_array[this[0]].spec1d_flux.value)
                var_array.append(skysubtract_array[this[0]].spec1d_uncertainty.quantity.value**2)
                mask_array.append(mask0)
                sky.append(sky_array[this[0]].spec1d_flux)
                if len(wavcal_array[this[0]].wav)>0:
                    wav_array.append(wavcal_array[this[0]].wav)
                else:
                    wav_array.append(np.full(len(skysubtract_array[0].spec1d_pixel),0.,dtype=float))
            else:
                data_array.append(np.full(len(skysubtract_array[0].spec1d_pixel),0.,dtype=float))
                var_array.append(np.full(len(skysubtract_array[0].spec1d_pixel),1.e+30,dtype=float))
                mask_array.append(np.full(len(skysubtract_array[0].spec1d_pixel),True,dtype=bool))
                sky.append(np.full(len(skysubtract_array[0].spec1d_pixel),0.,dtype=float))
                wav_array.append(np.full(len(skysubtract_array[0].spec1d_pixel),0.,dtype=float))

        var_array = np.array(var_array)
        mask_array = np.array(mask_array)
        data_array = np.array(data_array)
        sky = np.array(sky)
        wav_array = np.array(wav_array)

        ivar = 1./var_array
        mask_or = np.zeros((len(mask_array),len(mask_array[0])),dtype=int)
        for xxx in range(0,len(mask_or)):
            for yyy in range(0,len(mask_or[xxx])):
                if mask_array[xxx][yyy]:
                    mask_or[xxx][yyy] = 1
        mask_and = mask_or
                        
        if len(plugmap)!=len(skysubtract_array):
            print('error - plugmap has different length than data array!!!')
            print(len(plugmap),' vs',len(skysubtract_array))
            np.pause()
    else:
        wav_array = np.full(0,0.)
        data_array = np.full(0,0.)
        ivar = np.full(0,0.)
        mask_and = np.full(0,1)
        mask_or = np.full(0,1)
        sky = np.full(0,0.)

    primary_hdu = fits.PrimaryHDU(wav_array,header=data.header)
    new_hdul = fits.HDUList([primary_hdu])
    new_hdul.append(fits.ImageHDU(data_array,name='sky_subtracted'))
    new_hdul.append(fits.ImageHDU(ivar,name='inverse_variance'))
    new_hdul.append(fits.ImageHDU(mask_and,name='mask_and'))
    new_hdul.append(fits.ImageHDU(mask_or,name='mask_or'))
    new_hdul.append(fits.BinTableHDU(table_hdu,name='bin_table'))
    new_hdul.append(fits.ImageHDU(sky,name='mean_sky'))
    return new_hdul


def getfromfits(hdul):
    import numpy as np
    from astropy import time, coordinates as coord, units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    from astropy.io import fits

    class fitsobject:
        def __init__(self,radeg=None,decdeg=None,wav=None,spec=None,var=None,mask=None,sky_spec=None,obj=None,icode=None,
                     mjd=None,hjd=None,snratio=None,vheliocorr=None,header=None,filtername=None,channel=None,aperture=None,
                     run_id=None,field_name=None,wav_npoints=None,wav_rms=None,wav_resolution=None,wav_min=None,wav_max=None,
                     row=None,temperature=None):
            self.radeg = radeg
            self.decdeg = decdeg
            self.wav = wav
            self.spec = spec
            self.var = var
            self.mask = mask
            self.sky_spec = sky_spec
            self.obj = obj
            self.icode = icode
            self.mjd = mjd
            self.hjd = hjd
            self.snratio = snratio
            self.vheliocorr = vheliocorr
            self.header = header
            self.filtername = filtername
            self.channel = channel
            self.aperture = aperture
            self.run_id = run_id
            self.field_name = field_name
            self.wav_npoints = wav_npoints
            self.wav_rms = wav_rms
            self.wav_resolution=wav_resolution
            self.wav_min=wav_min
            self.wav_max=wav_max
            self.row=row
            self.temperature = temperature
#    nddata = astropy.nddata.CCDData.read(fitsfile[i],unit=u.electron)
#    hdul = fits.open(fitsfile)
    # for new CfA pipeline, HDUs are as follows (via email from Nelson on Apr 3, 2019:
    # hdu0: wavelength(angs)
    # hdu1: sky-subtracted, variance-weighted coadded spectra (total counts) OR flux-calibrated averaged spectra
    # hdu2: inverse variance (counts)
    # hdu3: AND bad pixel mask, requires pixel to be masked in all co-added frames
    # hdu4: OR bad pixel mask, requires pixel to be masked in any co-added frames
    # hdu5: plugmap structure (fiber info)
    # hdu6: combined sky spectra, absent if flux calibration was set
    # hdu7: summed (unweighted) spectra, absent if flux calibration was set
    wav = hdul[0].data
    skysub = hdul[1].data
    ivar = hdul[2].data
    mask_and = hdul[3].data
    mask_or = hdul[4].data
    mask = mask_or  # use the AND mask
    sky_spec = hdul[6].data
    var = 1./ivar
    masked = np.where(mask_or==1)[0]

    fiber = hdul[5].data
    obj = fiber['OBJTYPE']
    icode = fiber['ICODE']

    radeg = []
    decdeg = []
    mjd = []
    hjd = []
    vheliocorr = []
    snratio = []
    aperture = []
    channel = []
    filtername = []
    run_id = []
    field_name = []
    wav_npoints=[]
    wav_rms=[]
    wav_resolution=[]
    wav_min=[]
    wav_max=[]
    row=[]
    temperature=[]
    for j in range(0,len(obj)):
        ra,dec = fiber[j]['RA'],fiber[j]['DEC']
        radeg.append(ra)
        decdeg.append(dec)
        mjd.append(fiber[j]['MJD'])
        hjd.append(fiber[j]['HJD'])
        vheliocorr.append(fiber[j]['VHELIOCORR'])
        snratio.append(fiber[j]['SNRATIO'])
        aperture.append(fiber[j]['APERTURE'])
        channel.append(fiber[j]['CHANNEL'])
        filtername.append(fiber[j]['RESOLUTION'])
        run_id.append(fiber[j]['RUN_ID'])
        field_name.append(fiber[j]['FIELD_NAME'])
        wav_npoints.append(fiber[j]['WAV_NPOINTS'])
        wav_rms.append(fiber[j]['WAV_RMS'])
        wav_resolution.append(fiber[j]['WAV_RESOLUTION'])
        wav_min.append(fiber[j]['WAV_MIN'])
        wav_max.append(fiber[j]['WAV_MAX'])
        row.append(fiber[j]['ROW'])
        temperature.append(fiber[j]['TEMPERATURE'])
        if obj[j]=='TARGET':
            coords = coord.SkyCoord(ra,dec,unit=(u.deg,u.deg),frame='icrs')
            lco = coord.EarthLocation.from_geodetic(lon=-70.6919444*u.degree,lat=-29.0158333*u.degree,height=2380.*u.meter)
            times = time.Time(np.float(hdul[0].header['weighted_mjd']),format='mjd',location=lco)
            light_travel_time_helio = times.light_travel_time(coords,'heliocentric')
        bad = np.where(mask_or[j]==1)[0]
        var[j][bad] = 1.e+30
    radeg = np.array(radeg)
    decdeg = np.array(decdeg)
    mjd = np.array(mjd)
    hjd = np.array(hjd)
    snratio = np.array(snratio)
    vheliocorr = np.array(vheliocorr)
    aperture = np.array(aperture)
    channel = np.array(channel)
    filtername = np.array(filtername)
    run_id = np.array(run_id)
    field_name = np.array(field_name)
    wav_npoints = np.array(wav_npoints,dtype=np.object_)
    wav_rms = np.array(wav_rms,dtype=np.object_)
    wav_resolution = np.array(wav_resolution,dtype=np.object_)
    wav_min = np.array(wav_min,dtype=np.object_)
    wav_max = np.array(wav_max,dtype=np.object_)
    row = np.array(row)
    temperature = np.array(temperature,dtype=np.object_)
    header = hdul[0].header

    return fitsobject(radeg=radeg,decdeg=decdeg,wav=wav,spec=skysub,var=var,mask=mask,sky_spec=sky_spec,obj=obj,
                      icode=icode,mjd=mjd,hjd=hjd,snratio=snratio,vheliocorr=vheliocorr,header=header,filtername=filtername,
                      channel=channel,aperture=aperture,run_id=run_id,field_name=field_name,wav_npoints=wav_npoints,
                      wav_rms=wav_rms,wav_resolution=wav_resolution,wav_min=wav_min,wav_max=wav_max,row=row,temperature=temperature)


def multinest(fit_directory,root,targets,fitsobject,npix):
    import numpy as np
    from scipy.stats import skew,kurtosis
    from astropy.io import fits
    import os
    from os import path
    import glob

    class multinest_result:
        def __init__(self,posterior_1000=None,moments=None,bestfit_wav=None,bestfit_fit=None):
            self.posterior_1000 = posterior_1000
            self.moments = moments
            self.bestfit_wav = bestfit_wav
            self.bestfit_fit = bestfit_fit

    posterior_1000 = []
    moments = []
    bestfit_wav = []
    bestfit_fit = []

    for j in range(0,len(root)):
#    for j in targets:
        multinest_in = fit_directory+root[j]+'post_equal_weights.dat'
        bestfit_in = fit_directory+root[j]+'_bestfit.dat'

        bestfit_wav0 = np.zeros(npix)
        bestfit_fit0 = np.zeros(npix)
        bestfit_counts0 = np.zeros(npix)
        bestfit_varcounts0 = np.zeros(npix)

        posterior0_1000 = []
        for n in range(0,1000):
            posterior0_1000.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
#        posterior0_1000 = posterior0[np.random.randint(low=0,high=len(posterior0),size=1000)]
        moments0 = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

#        print(bestfit_in,path.exists(multinest_in),path.exists(bestfit_in))
        if ((path.exists(multinest_in))&(path.exists(bestfit_in))):
            with open(multinest_in) as f:
                data = f.readlines()
            posterior0 = []
            for line in data:
                p = line.split()
                posterior0.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),float(p[5]),
                                   float(p[6]),float(p[7]),float(p[8]),float(p[9]),float(p[10]),float(p[11]),
                                   float(p[12]),float(p[13]),float(p[14]),float(p[15])])
            posterior0 = np.array(posterior0)
            posterior0_1000 = posterior0[np.random.randint(low=0,high=len(posterior0),size=1000)]
#        moments0 = [[np.mean(posterior0[:,i]),np.std(posterior0[:,i]),skew(posterior0[:,i]),kurtosis(posterior0[:,i])] for i in range(0,len(p))]
            moments0 = [[np.mean(posterior0[:,i]) for i in range(0,len(p))],[np.std(posterior0[:,i]) for i in range(0,len(p))],[skew(posterior0[:,i]) for i in range(0,len(p))],[kurtosis(posterior0[:,i]) for i in range(0,len(p))]]

            with open(bestfit_in) as f:
                data = f.readlines()
            bestfit_wav0 = []
            bestfit_fit0 = []
            bestfit_counts0 = []
            bestfit_varcounts0 = []
            for line in data:
                p = line.split()
                bestfit_wav0.append(float(p[0]))
                bestfit_fit0.append(float(p[1]))
                bestfit_counts0.append(float(p[2]))
                bestfit_varcounts0.append(float(p[3]))
            bestfit_wav0 = np.array(bestfit_wav0)
            bestfit_fit0 = np.array(bestfit_fit0)
            bestfit_counts0 = np.array(bestfit_counts0)
            bestfit_varcounts0 = np.array(bestfit_varcounts0)

        bestfit_wav.append(bestfit_wav0)
        bestfit_fit.append(bestfit_fit0)

        posterior_1000.append(posterior0_1000)
        moments.append(moments0)

    posterior_1000 = np.array(posterior_1000)
    moments = np.array(moments)
    bestfit_wav = np.array(bestfit_wav)
    bestfit_fit = np.array(bestfit_fit)

#    if  len(bestfit_fit)!=len(targets):
#        print('ERROR: number of targets does not equal number of posterior files!!!!')
#        np.pause()

    return multinest_result(posterior_1000=posterior_1000,moments=moments,bestfit_wav=bestfit_wav,bestfit_fit=bestfit_fit)



def multinest_noflat(fit_directory,root,targets,fitsobject,npix):
    import numpy as np
    from scipy.stats import skew,kurtosis
    from astropy.io import fits
    import os
    from os import path
    import glob

    class multinest_result:
        def __init__(self,posterior_1000=None,moments=None,bestfit_wav=None,bestfit_fit=None):
            self.posterior_1000 = posterior_1000
            self.moments = moments
            self.bestfit_wav = bestfit_wav
            self.bestfit_fit = bestfit_fit

    posterior_1000 = []
    moments = []
    bestfit_wav = []
    bestfit_fit = []

    for j in range(0,len(root)):
#    for j in targets:
        multinest_in = fit_directory+root[j]+'post_equal_weights.dat'
        bestfit_in = fit_directory+root[j]+'_bestfit.dat'

        bestfit_wav0 = np.zeros(npix)
        bestfit_fit0 = np.zeros(npix)
        bestfit_counts0 = np.zeros(npix)
        bestfit_varcounts0 = np.zeros(npix)

        posterior0_1000 = []
        for n in range(0,1000):
            posterior0_1000.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
#        posterior0_1000 = posterior0[np.random.randint(low=0,high=len(posterior0),size=1000)]
        moments0 = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

#        print(bestfit_in,path.exists(multinest_in),path.exists(bestfit_in))
        if ((path.exists(multinest_in))&(path.exists(bestfit_in))):
            with open(multinest_in) as f:
                data = f.readlines()
            posterior0 = []
            for line in data:
                p = line.split()
                posterior0.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),float(p[5]),
                                   float(p[6]),float(p[7]),float(p[8]),float(p[9]),float(p[10]),float(p[11]),
                                   float(p[12]),float(p[13]),float(p[14]),float(p[15])])
            posterior0 = np.array(posterior0)
            posterior0_1000 = posterior0[np.random.randint(low=0,high=len(posterior0),size=1000)]
#        moments0 = [[np.mean(posterior0[:,i]),np.std(posterior0[:,i]),skew(posterior0[:,i]),kurtosis(posterior0[:,i])] for i in range(0,len(p))]
            moments0 = [[np.mean(posterior0[:,i]) for i in range(0,len(p))],[np.std(posterior0[:,i]) for i in range(0,len(p))],[skew(posterior0[:,i]) for i in range(0,len(p))],[kurtosis(posterior0[:,i]) for i in range(0,len(p))]]

            with open(bestfit_in) as f:
                data = f.readlines()
            bestfit_wav0 = []
            bestfit_fit0 = []
            bestfit_counts0 = []
            bestfit_varcounts0 = []
            for line in data:
                p = line.split()
                bestfit_wav0.append(float(p[0]))
                bestfit_fit0.append(float(p[1]))
                bestfit_counts0.append(float(p[2]))
                bestfit_varcounts0.append(float(p[3]))
            bestfit_wav0 = np.array(bestfit_wav0)
            bestfit_fit0 = np.array(bestfit_fit0)
            bestfit_counts0 = np.array(bestfit_counts0)
            bestfit_varcounts0 = np.array(bestfit_varcounts0)

#        else:
#            if not path.exists(multinest_in):
#                print(multinest_in)
#                print(bestfit_in)
#                cut1 = multinest_in.split('hjd')
#                look = cut1[0]+'hjd'+cut1[1][0:9]
#                os.system('ls '+look+'*post_equal_weights.dat')
#                exist = glob.glob(look+'*post_equal_weights.dat')
#                print(look+'*weights.dat')
#                print(exist)
#                np.pause()


        bestfit_wav.append(bestfit_wav0)
        bestfit_fit.append(bestfit_fit0)

        posterior_1000.append(posterior0_1000)
        moments.append(moments0)

    posterior_1000 = np.array(posterior_1000)
    moments = np.array(moments)
    bestfit_wav = np.array(bestfit_wav)
    bestfit_fit = np.array(bestfit_fit)

#    if  len(bestfit_fit)!=len(targets):
#        print('ERROR: number of targets does not equal number of posterior files!!!!')
#        np.pause()

    return multinest_result(posterior_1000=posterior_1000,moments=moments,bestfit_wav=bestfit_wav,bestfit_fit=bestfit_fit)



def get_plot_resolution(thar,root0):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
#    matplotlib.use('TkAgg')
    matplotlib.use('pdf')

    plot_resolution_file = root0+'_resolution.pdf'
    gs = plt.GridSpec(7,7)   #  define multi-panel plot
    gs.update(wspace=0,hspace=0)   #  specify inter-panel spacing
    fig = plt.figure(figsize=(6,6))   #  define plot size
    ax1 = fig.add_subplot(gs[0:6,0:6])

    for j in range(0,len(thar)):
        for k in range(0,len(thar[j])):
            if len(thar[j][k].wav)>1:
                keep = np.where(thar[j][k].wav.mask==False)[0]
                if len(keep)>1:
                    dwav = np.max(thar[j][k].wav[keep])-np.min(thar[j][k].wav[keep])
                    dpix = np.max([thar[j][k].fit_lines.fit[q].mean.value for q in keep])-np.min([thar[j][k].fit_lines.fit[q].mean.value for q in keep])
#                                ax1.plot([thar[j][k].fit_lines.fit[q].mean.value for q in keep],thar[j][k].resolution([thar[j][k].fit_lines.fit[q].mean.value for q in keep])*dwav/dpix,alpha=0.2,lw=0.3)
                    ax1.plot([thar[j][k].wav[q] for q in keep],thar[j][k].resolution([thar[j][k].fit_lines.fit[q].mean.value for q in keep])*dwav/dpix,alpha=0.2,lw=0.3)
                                                    

    ax1.set_xlabel('$\lambda$ [Angs.]')
    ax1.set_ylabel('$\sigma_{\lambda}$ [Angs.]')
    plt.savefig(plot_resolution_file,dpi=200)
#    plt.show()
    plt.close()
    return


def get_mapflat_compare(apertures_profile_middle,map_apertures_profile_middle,mapflat_compare_iterations,mapflat_compare_order):
    import numpy as np
    import scipy
    from astropy.modeling import models,fitting
    import matplotlib.pyplot as plt

    keep = np.where(apertures_profile_middle.realvirtual)[0]
    x = np.array([apertures_profile_middle.fit[q].mean.quantity.value for q in keep])
    y = np.ma.masked_array([map_apertures_profile_middle.fit[q].mean.quantity.value for q in keep],mask=np.full(len(x),False,dtype='bool'))
    bad = np.where(map_apertures_profile_middle.realvirtual==False)[0]
    y.mask[bad] = True
    map_flat_compare_iterations = 10
    mapflat_compare_order = 4
    mapflat_compare_init = models.Polynomial1D(degree=mapflat_compare_order)
    mapflat_compare_fitter = fitting.LinearLSQFitter()
    mapflat_compare_func = mapflat_compare_init
    use = np.where(y.mask==False)[0]
    mapflat_compare_rms = np.sqrt(np.mean((y[use]-mapflat_compare_func(x[use]))**2))
    for q in range(0,map_flat_compare_iterations):
        use = np.where(y.mask==False)[0]
        if len(use)>3.:
            mapflat_compare_func = mapflat_compare_fitter(mapflat_compare_init,x[use],y[use])
            mapflat_compare_rms = np.sqrt(np.mean((y[use]-mapflat_compare_func(x[use]))**2))
            outlier = (np.where(np.abs(y-mapflat_compare_func(x))>3.*mapflat_compare_rms))[0]    # reject outliers
            y.mask[outlier] = True
    use = np.where(y.mask==False)[0]
    mapflat_compare_npoints = len(use)
    return mapflat_compare(mapflat_compare_func,mapflat_compare_rms,mapflat_compare_npoints)

