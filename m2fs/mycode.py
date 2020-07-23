# etaxi module
import numpy as np
import string

def weightedmean(x,sigx):
    sum1=np.sum(x/sigx**2)
    sum2=np.sum(1./sigx**2)
    return sum1/sum2,np.sqrt(1./sum2)

def meansigma2_maxl(x,ex,mean0,sigma20,iterations):
    for i in range(0,iterations):
        sum1=0.
        sum2=0.
        sum3=0.
        sum4=0.
        for j in range(0,len(x)):
            sum1=sum1+x[j]/(1.+ex[j]**2/sigma20**2)
            sum2=sum2+1./(1.+ex[j]**2/sigma20**2)
            sum3=sum3+(x[j]-mean0)**2/(1.+ex[j]**2/sigma20**2)**2
        mean1=sum1/sum2
        sigma21=sum3/sum2
        mean0=mean1
        sigma20=sigma21
    return mean1,sigma21

def etaxiarr(rarad,decrad,racenter,deccenter):
    xi=np.zeros(len(rarad))
    eta=np.zeros(len(rarad))
    for i in range(len(rarad)):
        xi[i]=np.cos(decrad[i])*np.sin(rarad[i]-racenter[i])/(np.sin(deccenter[i])*np.sin(decrad[i])+np.cos(deccenter[i])*np.cos(decrad[i])*np.cos(rarad[i]-racenter[i]))*180.*60./np.pi
        eta[i]=(np.cos(deccenter[i])*np.sin(decrad[i])-np.sin(deccenter[i])*np.cos(decrad[i])*np.cos(rarad[i]-racenter[i]))/(np.sin(deccenter[i])*np.sin(decrad[i])+np.cos(deccenter[i])*np.cos(decrad[i])*np.cos(rarad[i]-racenter[i]))*180.*60./np.pi
#        print rarad[i],decrad[i],xi[i],eta[i]
    return xi,eta

def radecradians(rah,ram,ras,chardecd,decm,decs):
    rarad=(rah+ram/60.+ras/3600.)*360./24.*np.pi/180.
    [signdec,decd]=decsign(chardecd)
    decrad=(float(decd)+float(decm)/60.+decs/3600.)*np.pi/180.
    if signdec == '-' :
        decrad=-1.*decrad
    return rarad,decrad

def radecradiansarr(rah,ram,ras,chardecd,decm,decs):
    rarad=np.zeros(len(ras))
    decrad=np.zeros(len(ras))
    for i in range(len(ras)):
        rarad[i]=(rah[i]+ram[i]/60.+ras[i]/3600.)*360./24.*np.pi/180.
        [signdec,decd]=decsign(chardecd[i])
        decrad[i]=(float(decd)+float(decm[i])/60.+decs[i]/3600.)*np.pi/180.
        if signdec == '-' :
            decrad[i]=-1.*decrad[i]
    return rarad,decrad

def decsign(chardecd):
    if '-' in chardecd:
        signdec='-'
        decd=int((string.split(chardecd,'-'))[1])
    else:
        signdec='+'
        decd=int(chardecd)
    return signdec,decd

def paradarr(x,y):
    pa=np.zeros(len(x))
    for i in range(len(x)):
        if x[i] >= 0.:
            if y[i] >= 0.:
                pa[i]=(np.pi/2.-np.arctan(y[i]/x[i]))*180./np.pi
            if y[i] < 0.:
                pa[i]=(np.pi-np.arctan(-1.*x[i]/y[i]))*180./np.pi
        if x[i] < 0.:
            if y[i] >= 0.:
                pa[i]=(3.*np.pi/2.+np.arctan(-1.*y[i]/x[i]))*180./np.pi
            if y[i] < 0.:
                pa[i]=(np.pi+np.arctan(x[i]/y[i]))*180./np.pi
    return pa

def calc_min_interval(x, alpha):
    """Internal method to determine the minimum interval of
    a given width
    Assumes that x is sorted numpy array.
    """

    n = len(x)
    cred_mass = 1.0 - alpha

    interval_idx_inc = int(np.floor(cred_mass * n))
    n_intervals = n - interval_idx_inc
    interval_width = x[interval_idx_inc:] - x[:n_intervals]

    if len(interval_width) == 0:
        print_('Too few elements for interval calculation')
        return None, None

    min_idx = np.argmin(interval_width)
    hdi_min = x[min_idx]
    hdi_max = x[min_idx + interval_idx_inc]
    return hdi_min, hdi_max

def hpd(x, alpha):
    """Calculate highest posterior density (HPD) of array for given alpha. The HPD is the
    minimum width Bayesian credible interval (BCI).
    :Arguments:
      x : Numpy array
          An array containing MCMC samples
      alpha : float
          Desired probability of type I error
    """

    # Make a copy of trace
    x = x.copy()

    # For multivariate node
    if x.ndim > 1:

        # Transpose first, then sort
        tx = tr(x, list(range(x.ndim)[1:]) + [0])
        dims = shape(tx)

        # Container list for intervals
        intervals = np.resize(0.0, (2,) + dims[:-1])

        for index in make_indices(dims[:-1]):

            try:
                index = tuple(index)
            except TypeError:
                pass

            # Sort trace
            sx = np.sort(tx[index])

            # Append to list
            intervals[0][index], intervals[1][index] = calc_min_interval(sx, alpha)
        
        # Transpose back before returning
        return np.array(intervals)

    else:
        # Sort univariate node
        sx = np.sort(x)

        return np.array(calc_min_interval(sx, alpha))


def gal_uvw(distance=None, lsr=None, ra=None, dec=None, pmra=None, pmdec=None, vrad=None, plx=None):
   """
    NAME:
        GAL_UVW
    PURPOSE:
        Calculate the Galactic space velocity (U,V,W) of star
    EXPLANATION:
        Calculates the Galactic space velocity U, V, W of star given its
        (1) coordinates, (2) proper motion, (3) distance (or parallax), and
        (4) radial velocity.
    CALLING SEQUENCE:
        GAL_UVW [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD= , DISTANCE=
                 PLX= ]
    OUTPUT PARAMETERS:
         U - Velocity (km/s) positive toward the Galactic *anti*center
         V - Velocity (km/s) positive in the direction of Galactic rotation
         W - Velocity (km/s) positive toward the North Galactic Pole
    REQUIRED INPUT KEYWORDS:
         User must supply a position, proper motion,radial velocity and distance
         (or parallax).    Either scalars or vectors can be supplied.
        (1) Position:
         RA - Right Ascension in *Degrees*
         Dec - Declination in *Degrees*
        (2) Proper Motion
         PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
         PMDEC = Proper motion in Declination (typically mas/yr)
        (3) Radial Velocity
         VRAD = radial velocity in km/s
        (4) Distance or Parallax
         DISTANCE - distance in parsecs
                    or
         PLX - parallax with same distance units as proper motion measurements
               typically milliarcseconds (mas)
   
    OPTIONAL INPUT KEYWORD:
         /LSR - If this keyword is set, then the output velocities will be
                corrected for the solar motion (U,V,W)_Sun = (-8.5, 13.38, 6.49)
                (Coskunoglu et al. 2011 MNRAS) to the local standard of rest.
                Note that the value of the solar motion through the LSR remains
                poorly determined.
     EXAMPLE:
         (1) Compute the U,V,W coordinates for the halo star HD 6755.
             Use values from Hipparcos catalog, and correct to the LSR
         ra = ten(1,9,42.3)*15.    & dec = ten(61,32,49.5)
         pmra = 627.89  &  pmdec = 77.84         ;mas/yr
         dis = 144    &  vrad = -321.4
         gal_uvw,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,dis=dis,/lsr
             ===>  u=154  v = -493  w = 97        ;km/s
   
         (2) Use the Hipparcos Input and Output Catalog IDL databases (see
         http://idlastro.gsfc.nasa.gov/ftp/zdbase/) to obtain space velocities
         for all stars within 10 pc with radial velocities > 10 km/s
   
         dbopen,'hipparcos,hic'      ;Need Hipparcos output and input catalogs
         list = dbfind('plx>100,vrad>10')      ;Plx > 100 mas, Vrad > 10 km/s
         dbext,list,'pmra,pmdec,vrad,ra,dec,plx',pmra,pmdec,vrad,ra,dec,plx
         ra = ra*15.                 ;Need right ascension in degrees
         GAL_UVW,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,plx = plx
         forprint,u,v,w              ;Display results
    METHOD:
         Follows the general outline of Johnson & Soderblom (1987, AJ, 93,864)
         except that U is positive outward toward the Galactic *anti*center, and
         the J2000 transformation matrix to Galactic coordinates is taken from
         the introduction to the Hipparcos catalog.
    REVISION HISTORY:
         Written, W. Landsman                       December   2000
         fix the bug occuring if the input arrays are longer than 32767
           and update the Sun velocity           Sergey Koposov June 2008
   	   vectorization of the loop -- performance on large arrays
           is now 10 times higher                Sergey Koposov December 2008
   """

   n_params = 3
   
   if n_params == 0:   
      print( 'Syntax - GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD=')
      print( '                  Distance=, PLX=')
      print( '         U, V, W - output Galactic space velocities (km/s)')
      return None
   
   if ra is None or dec is None:   
      raise Exception('ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)')
   if plx is None and distance is None:
      raise Exception('ERROR - Either a parallax or distance must be specified')
   if distance is not None:
      if np.any(distance==0):
         raise Exception('ERROR - All distances must be > 0')
      plx = 1e3 / distance          #Parallax in milli-arcseconds
   if plx is not None and np.any(plx==0):   
      raise Exception('ERROR - Parallaxes must be > 0')
   
   cosd = np.cos(np.deg2rad(dec))
   sind = np.sin(np.deg2rad(dec))
   cosa = np.cos(np.deg2rad(ra))
   sina = np.sin(np.deg2rad(ra))
   
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s   
   a_g = np.array([[0.0548755604, +0.4941094279, -0.8676661490],
                [0.8734370902, -0.4448296300, -0.1980763734], 
                [0.4838350155, 0.7469822445, +0.4559837762]])
   
   vec1 = vrad
   vec2 = k * pmra / plx
   vec3 = k * pmdec / plx
   
   u = (a_g[0,0] * cosa * cosd + a_g[1,0] * sina * cosd + a_g[2,0] * sind) * vec1 + (-a_g[0,0] * sina + a_g[1,0] * cosa) * vec2 + (-a_g[0,0] * cosa * sind - a_g[1,0] * sina * sind + a_g[2,0] * cosd) * vec3
   v = (a_g[0,1] * cosa * cosd + a_g[1,1] * sina * cosd + a_g[2,1] * sind) * vec1 + (-a_g[0,1] * sina + a_g[1,1] * cosa) * vec2 + (-a_g[0,1] * cosa * sind - a_g[1,1] * sina * sind + a_g[2,1] * cosd) * vec3
   w = (a_g[0,2] * cosa * cosd + a_g[1,2] * sina * cosd + a_g[2,2] * sind) * vec1 + (-a_g[0,2] * sina + a_g[1,2] * cosa) * vec2 + (-a_g[0,2] * cosa * sind - a_g[1,2] * sina * sind + a_g[2,2] * cosd) * vec3

   lsr_vel = np.array([-8.5, 13.38, 6.49]) # notice the sign of the first velocity 
   											#component, it is negative because 
									# in this program U points toward anticenter
#   if (lsr is not None):   
#       print 'bbbbbb'
   u = u + lsr_vel[0]
   v = v + lsr_vel[1]
   w = w + lsr_vel[2]
   
   return u,v,w

def oldstats(v,sigv,p):
    count=np.size(v)
    vdisp=100.
    if(count==1):
        vmean=v[0]
        vdisp=-999.
        sigvmean=-999.
        sigvdisp=999.
        return mean,vdisp,sigvmean,sigvdisp
    vmean=np.mean(v)
    for i in range(0,50):
        sum1=np.sum(p*v/(1.+sigv**2/vdisp**2))
        sumweight=np.sum(p/(1.+sigv**2/vdisp**2))
        vmean=sum1/sumweight
        sum2=np.sum(p*((v-vmean)**2)/(1.+sigv**2/vdisp**2)**2)
        vdisp=np.sqrt(sum2/sumweight)
        d1=np.sum(-p/(vdisp**2+sigv**2))
        d2=np.sum(-p/(vdisp**2+sigv**2)+2.*vdisp**2*p/(vdisp**2+sigv**2)**2+p*(v-vmean)**2/(vdisp**2+sigv**2)**2-4.*vdisp**2*p*(v-vmean)**2/(vdisp**2+sigv**2)**3)
        d3=np.sum(-2.*p*vdisp**2*(v-vmean)/(sigv**2+vdisp**2)**2)
        a=d2/(d1*d2-d3**2)
        b=d1/(d1*d2-d3**2)
        sigvmean=np.sqrt(abs(a))
        sigvdisp=np.sqrt(abs(b))
    return vmean,vdisp,sigvmean,sigvdisp

def stdmean(a,axis=None,dtype=None,out=None,ddof=0,keepdims=np._NoValue):
    std=np.ma.std(a)
#    print(len(a))
    return std/np.sqrt(len(a))
