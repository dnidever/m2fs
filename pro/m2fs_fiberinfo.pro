pro m2fs_fiberinfo,head,fiberstr

; Extract M2Fs fiber info from header

; Not enough inputs
if n_elements(head) eq 0 then begin
  print,'Syntax - m2fs_fiberinfo,head,fiberstr'
  return
endif

gd = where(strmid(head,0,5) eq 'FIBER',ngd)
print,strtrim(ngd,2),' fibers found in header'

fiberstr = replicate({headid:0L,type:'',sra:'',sdec:'',ra:0.0d,dec:0.0d0,fid:0L,fab:'',ch:'',line:''},ngd)

for i=0,ngd-1 do begin
  line = head[gd[i]]

  ;FIBER101= '14:00:46.52_-42:44:25.0'    / FID=125 Fab=24/1 C/H=02/01   

  fiberstr[i].line = line
  headid = strmid(line,5,3)
  fiberstr[i].headid = long(headid)
  lo = strpos(line,'=')
  hi = strpos(line,'/')
  val = strtrim(strmid(line,lo+1,hi-lo-1),2)
  ; Remove beg/end tick marks if they are there
  if strmid(val,0,1) eq "'" then val=strmid(val,1)
  if strmid(val,strlen(val)-1,1) eq "'" then val=strmid(val,0,strlen(val)-1)
  case val of
    'sky': begin
        fiberstr[i].type='sky'
      end    
    'inactive': begin
        fiberstr[i].type='inactive'
      end
    else: begin
        fiberstr[i].type='object'
        dum = strsplit(val,'_',/extract)
        fiberstr[i].sra = dum[0]
        fiberstr[i].sdec = dum[1]
        fiberstr[i].ra = sexig2ten(fiberstr[i].sra)*15.
        fiberstr[i].dec = sexig2ten(fiberstr[i].sdec)
      end
  endcase

  ; FID, Fab and C/H values in the comment
  com = strtrim(strmid(line,hi+1),2)
  dum = strsplit(com,' ',/extract)
  fid = (strsplit(dum[0],'=',/extract))[1]
  fiberstr[i].fid = long(fid)
  fab = (strsplit(dum[1],'=',/extract))[1]
  fiberstr[i].fab = fab
  ch = (strsplit(dum[2],'=',/extract))[1]
  fiberstr[i].ch = ch
endfor

;stop

end
