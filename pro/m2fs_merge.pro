pro m2fs_merge,input,clobber=clobber

; Merge the four M2FS amp files into one
; perform overscan subtraction and trimming

; load input
loadinput,input,files,count=nfiles

dir = file_dirname(files)
base = file_basename(files)
prefix = strmid(base,0,5)  ; r/b and 4 digit number

; Get unique names
;  allow for the possibility that there are files
;  in separate directories but with the same name
temp = dir+'/'+prefix
ui = uniq(temp,sort(temp))
udir = dir[ui]
uprefix = prefix[ui]
nframes = n_elements(ui)

; Loop through the frames
for i=0LL,nframes-1 do begin

  print,strtrim(i+1,2),'/',strtrim(nframes,2),' ',uprefix[i]
   
  ampfiles = udir[i]+'/'+uprefix[i]+'c'+strtrim(indgen(4)+1,2)+'.fits'
  test = file_test(ampfiles)
  if total(test) lt 4 then begin
    bd = where(test eq 0,nbd)
    print,ampfiles[bd],' NOT FOUND'
    goto,BOMB
  endif

  ; Check if outfile exists
  outfile = udir[i]+'/'+uprefix[i]+'.fits'
  if file_test(outfile) and not keyword_set(clobber) then begin
    print,'Outfile ',outfile,' already exists and /clobber not set.  Skipping.'
    goto,BOMB 
  endif
  
  ;Mario has written a short set of commands to pack separate amplifer images 
  ;together.  Put these in a file, one set for EVERY SINGLE IMAGE and feed that
  ;to the IRAF command line (sample contents of file, packr.cl, below):
  ;
  ;imcreate r1010.fits 2 4096 4112 header=copy pixtype=real reference='r1010c1.fits'
  ;imcopy r1010c1.fits[1:2048,1:2056] r1010.fits[1:2048,4112:2057]
  ;imcopy r1010c2.fits[1:2048,1:2056] r1010.fits[4096:2049,4112:2057]
  ;imcopy r1010c3.fits[1:2048,1:2056] r1010.fits[4096:2049,1:2056]
  ;imcopy r1010c4.fits[1:2048,1:2056] r1010.fits[1:2048,1:2056]

  ; Loop through the chips
  for j=0,3 do begin
  
    fits_read,ampfiles[j],im,head,/no_abort,message=message
    im = float(im)

    if j eq 0 then begin
      binning = sxpar(head,'binning',count=nbinning)
      if nbinning gt 0 then begin
        dum = strsplit(binning,'x',/extract)
        xbin = long(dum[0])
        ybin = long(dum[1])
        print,'Binning = ',binning
      endif else begin
        print,'No binning'
        xbin = 1
        ybin = 1
      endelse
      nx = 2048/xbin
      ny = 2048/ybin
      nxf = 4096/xbin               ; final size
      nyf = 4096/ybin
      fim = fltarr(nxf,nyf)
      fhead = head  
      sxaddpar,fhead,'bitpix',-32
      sxaddpar,fhead,'naxis1',nxf
      sxaddpar,fhead,'naxis1',nyf
      sxaddpar,fhead,'bscale',1.0
      sxaddpar,fhead,'bzero',0.0
    endif
      
    ; Can use BIASSEC, DATASEC, TRIMSEC
  
    ; Overscan is on the right-hand side and the top
    oim = im[2048/xbin:*,0:2048/ybin-1]
    omed = median(oim,dim=1)
    omedbin = rebin(omed,n_elements(omed)/64)
    xomedbin = lindgen(n_elements(omedbin))*64+64/2
    x = findgen(nx)
    interp,xomedbin,omedbin,x,oval
    ; maybe use B-splines instead
    sxaddhist,ampfiles[j]+' median overscan '+strtrim(median(omed),2),fhead
    print,'Median overscan ',strtrim(median(omed),2),' counts'

    ;plot,omed,/ysty
    ;oplot,xomedbin,omedbin,ps=-1,sym=2,co=250
    ;oplot,x,oval,co=200
    ;stop
    
    ; Overscan correct/subtract
    subim = im[0:nx-1,0:ny-1] - replicate(1,nx)#oval

    ; Add to final image
    ;imcopy r1010c1.fits[1:2048,1:2056] r1010.fits[1:2048,4112:2057]
    ;imcopy r1010c2.fits[1:2048,1:2056] r1010.fits[4096:2049,4112:2057]
    ;imcopy r1010c3.fits[1:2048,1:2056] r1010.fits[4096:2049,1:2056]
    ;imcopy r1010c4.fits[1:2048,1:2056] r1010.fits[1:2048,1:2056]
    case j of
    0: fim[0:nx-1,ny:2*ny-1] = reverse(subim,2)
    1: fim[nx:2*nx-1,ny:2*ny-1] = reverse(reverse(subim,1),2)
    2: fim[nx:2*nx-1,0:ny-1] = reverse(subim,1)
    3: fim[0:nx-1,0:ny-1] = subim
    else: stop
    endcase
  endfor

  ; Write to BASE.fits, i.e. r1234.fits
  print,'Writing merged file to ',outfile
  MWRFITS,fim,outfile,fhead,/create
  
  ;stop

  BOMB:

endfor


;stop

end
