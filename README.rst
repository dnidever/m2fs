M2FS Reduction Software
=======================

m2fs_zero.py
m2fs_dark.py
m2fs_reduce.py
m2fs_apertures.py

Run the above in that sequence.  All of them refer to m2fs_process.py

m2fs_zero
---------
Creates master bias frames.

m2fs_dark
---------
Creates master dark frames.

m2fs_reduce
-----------
Performs overscan/bias/dark corrections, stitches together individual amps into single frame with suffix '_stitched.fits'


m2fs_apertures
--------------

- Can always run many steps at the same time.
- The ``overwrite=True`` parameter sets whether any step will overwrite any existing data that it previously created.

Run the code like this: ``python m2fs_apertures.py RUN_science_raw``
For example: ``python m2fs_apertures.py nov2013_science_raw``

where RUN_science_raw is the input data for the fields:

ut20131123 0233 0241 0241 0239-0240-0236 Outer_LMC_Set3 0234-0235-0236-0237-0238 /net/dl1/users/dnidever/net/halo/dln5q/doradus/research/m2fs/data/nov2013/fibermap/outer_lmc_1.plate none lmc3


**RUN_science_raw:**

- column 1: UT date and directory name
- column 2: beginning frame number of the sequence
- column 3: end frame number of the sequence
- column 4: frame number of continuum-source calibration (eg. quartz lamp, LED, etc.)  
- column 5: frame number(s) of arc-lamp spectra (if more than one, delimit using hyphen (e.g., 0435-0439-0441)
- column 6: character string naming the field (I generally use what was written in observing log)
- column 7: frame number(s) of science exposures (if more than one, delimit using hyphen)
- column 8: location of .fibermap file that contains fiber allocations
- column 9: string that contains instructions for updating fiber allocations to account for any plugging errors documented in observing log (hopefully there are none, in which case write 'none')
- column 10: shorthand name of object (eg for any LMC field this becomes 'lmc')

**1) Trim**

- set trim=True, everything else to False
- set lower and upper boundaries interactively on quartz
- creates a _image_boundary.pickle file for each exposure
  
**2) Initialize**
   
- set initialize=True, everything else to False
- automatically finds the peaks, takes a while
- fits gaussian profile for each peak
- creates a _columnspec_array.pickle file for each exposure

**3) Find**
   
- set find=True, everything else to False
- interactive
- update/add/delete apertures
- assign artificial aperture that wasn't plugged
- "n" for new real aperture
- "a" for new artificial aperture
- all fibers are plugged in the twilight exposure
- need to deal with broken fibers as well
- need to identify 128 apertures
- creates a _apertures_profile_middle.pickle file
  
**4) Trace**
   
- set trace_all=True, everything else to False
- automated, takes some time as well
- uses data from "initialize" set to fit traces for each aperture
- some predefined parameters, should work okay
- can use trace_edit=True to manually edit, but generally not necessary
- there is code plot_trace() in m2fs_process.py
- uncomment code in m2fs_process.py lines 932-945 to see/interact with traces, hit 'q' to go to next one
- set make_image=True to get PSF of data+traces
- creates a _aperture_array.pickle and _apertures2d.pdf if you set make_image=True
  
**5) apmask**
   
- mask apertures
- set apmask=True, everything else to False
- automated, fairly fast
- creates _apmask.pickle file

**6) apflat**
   
- flat fielding
- set apflat=True, everything else to False
- automated
- using quart spectra for flat fields
- created _apflat.pickle and _apflat_residual.pickle files
  
**7) flat field correction**
   
- apply flat field correction
- set apflatcorr=True
- automated
- can run at same tiem as apflat
- creates _apflatcorr.pickle files
- can inspect output pickle files

**8) scattered light**
   
- set scatteredlightcorr=True
- fit 2D model to inter-aperture region and subtract it
- automated, fairly fast
- creates _scatteredlightcorr.pickle files
  
**9) extraction**
   
- done separately for flat, thar, sci
- they all use the apertures/traces from the quartz
- automated
- can set all three to run at same time
- creates _extract1d_array.pickle files

**10) ID lines template**
    
- generate template
- set id_lines_template=True
- similar to IRAF identify
- don't have a reference solution for my setup, I'll need to do it manually
- once I do this once I can use the template for other frames that use the same setup
- arc-specific template, one for ThNe, one for ThAr, and one for the skylines
- output is a pickle file that contains the wavelength solutions (wave vs. pixel) an function to fit that.  will have file number in it. rename it to something meaningfule, e.g. NeAr_idlines, ...
- creates _id_lines_template.pickle files
- only need to do this ONCE for a given setup/arc type.  once you have it you can copy the _id_lines_template.pickle file to the exposure name of the other arcs of the same type.
  
**11) ID lines translate**
    
- applies shifts and stretches (polynomial stretch function)
- give it template name
- looks for the information in "arc_templates"
- similar to IRAF reidentify
- fits template spectrum to arc spectrum
- shift and stretch fit
- automated
- uses "dynesty" package
- 1) uses dynesty to find the global minimum, like MCMC
- 2) gradient descent
- can fail a few % of the time
- likely an issue with tracing/extraction
- creates _id_lines_array.pickle files
  
**12) wavecal**

- set wavcal=True
- this performs the wavelength calibration of the science exposures using the wavelength solution of the arcs
- creates _wavcal_array.pickle files

**13) cr_reject**

- set cr_reject=True
- this performs cosmic ray rejection on the object exposures
- BEWARE: this can have problems with sky lines
- creates _cr_reject_array.pickle
- you can skip this step by copying the _extract1d_array.pickle files to _cr_reject_array.pickle
  
**13) stack_twilight**

- set stack_twilight=True
- only to be run on a sequence of twilight exposures
- this will stack the twilight subexposures
- used in the throughputcorr step for throughput correction of the object exposures
- creates _twilightstack_array.pickle and twilightstack_wavcal_array.pickle files with the name of the twilight, run, and spectrograph channel (e.g., b_twilight1_nov2013_twilightstack_wavcal_array.pickle, b_twilight1_nov2013_twilightstack_array.pickle)
  
**14) throughputcorr**

- set throughputcorr=True
- uses the stacked twilight frames to determine and apply fiber throughput corrections (as function of wavelength)
- creates _throughput_array.pickle and _throughputcorr_array.pickle files

**15) plugmap**

- set plugmap=True
- gets object information from the plugmap
- creats _plugmap.pickle files

**16) skysubtract**

- set skysubtract=True
- (not interactive) uses fibermap information to identify sky spectra, combines individual sky spectra to obtain mean sky spectrum and subtracts mean sky spectrum from individual spectra
- create _sky_array.pickle and _skysubtract_array.pickle files

**17) stack_frames**

- set stack_frames=True
- (not interactive) stacks science subexposures
- creates _stack_wavcal_array.pickle and _stack_array.pickle files with the name of the field, run, and spectrograph channel (e.g., b_Outer_LMC_Set3_nov2013_stack_wavcal_array.pickle and b_Outer_LMC_Set3_nov2013_stack_array.pickle)

**18) writefits**

- set writefits=True
- (not interactive) writes stacked science frames to fits files
- creates _skysubtract.fits files for each science frame (with a datestamp in the name, e.g. b0234_2013-11-23_04:58:31_skysubtract.fits), and also _stackskysub_file for the stacked frames with field name, run, spectrograph channel, and timestamp in the name (e.g., b_Outer_LMC_Set1_nov2013_2013-11-23_07:42:45_stackskysub_file).

