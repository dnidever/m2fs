M2FS Reduction Software
=======================

m2fs_zero.py
m2fs_dark.py
m2fs_reduce.py
m2fs_apertures.py

Run the above in that sequence.  All of them refer to m2fs_process.py

Procedures
----------

Can always run many steps at the same time.

**1) Trim**

- set trim=True, everything else to False
- set lower and upper boundaries interactively on quartz

**2) Initialize**
   
- set initialize=True, everything else to False
- automatically finds the peaks, takes a while
- fits gaussian profile for each peak

3) Find
   
- set find=True, everything else to False
- interactive
- update/add/delete apertures
- assign artificial aperture that wasn't plugged
- "n" for new real aperture
- "a" for new artificial aperture
- all fibers are plugged in the twilight exposure
- need to deal with broken fibers as well
- need to identify 128 apertures

4) Trace
   
- set trace_all=True, everything else to False
- automated, takes some time as well
- uses data from "initialize" set to fit traces for each aperture
- some predefined parameters, should work okay
- can use trace_edit=True to manually edit, but generally not necessary
- there is code plot_trace() in m2fs_process.py
- uncomment code in m2fs_process.py lines 932-945 to see/interact with traces, hit 'q' to go to next one
- set make_image=True to get PSF of data+traces

5) apmask
   
- mask apertures
- set apmask=True, everything else to False
- automated, fairly fast

6) apflat
   
- flat fielding
- set apflat=True, everything else to False
- automated
- using quart spectra for flat fields

7) flat field correction
   
- apply flat field correction
- set apflatcorr=True
- automated
- can run at same tiem as apflat
- can inspect output pickle files

8) scattered light
   
- set scatteredlightcorr=True
- fit 2D model to inter-aperture region and subtract it
- automated, fairly fast

9) extraction
   
- done separately for flat, thar, sci
- they all use the apertures/traces from the quartz
- automated
- can set all three to run at same time

10) ID lines template
    
- generate template
- set id_lines_template=True
- similar to IRAF identify
- don't have a reference solution for my setup, I'll need to do it manually
- once I do this once I can use the template for other frames that use the same setup
- arc-specific template, one for ThNe, one for ThAr, and one for the skylines
- output is a pickle file that contains the wavelength solutions (wave vs. pixel) an function to fit that.  will have file number in it. rename it to something meaningfule, e.g. NeAr_idlines, ...

11) ID lines translate
    
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





