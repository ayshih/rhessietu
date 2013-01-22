; NAME: ingest
;
; PURPOSE: Processes an SSR file from the RHESSI ETU/GSE
;
; INPUTS:
;   filename	SSR file from GSE
;   source	optional - specifies a source and passes it to "autoroi" (x2)
;   slice	optional - divides the spectrum into slices of this many seconds
;   plot	keyword - diagnostic plot of spectra and lightcurves
;   save	keyword - creates a save file, XXXX.ssr => XXXX.sav
;
; OUTPUTS:
;   spec_front	spectrum from the front segment (2D if slice is set)
;   spec_rear	spectrum from the rear segment (2D if slice is set)
;   time	total acquisition time (may be inaccurate from bad packets)
;   gain_front	optional - energy calibration for front if source is set
;   gain_rear	optional - energy calibration for rear if source is set
;
; EXAMPLES:
;   ingest,'12339202.ssr',sf,sr,time,/plot
;   ingest,'12339203.ssr',source=2
;
; HISTORY:
;   2012-12-09, AYS: release
;   2012-12-12, AYS: fix to rear segment source specification for autoroi
;   2012-12-14, AYS: only use the first filename supplied (or the first of a wildcard match)
;   2012-12-17, AYS: added plot indication of slices and optional outputs for gain
;   2013-01-21, AYS: removes duplicate reset events (using a window of 40 us)

pro ingest,filename,spec_front,spec_rear,time,$
  source=source,gain_front=gain_front,gain_rear=gain_rear,$
  slice=slice,plot=plot,save=save

ssrfile = (file_search(filename, count=count))[0]
if count eq 0 then begin
  print,"No such file"
  return
endif

o=hsi_eventlist()
d=o->getdata(filename=ssrfile,file_type='gse',time_range=[0,0])

; Known quirks with SSR GSE files
;   The timer rolls over after 2^32 binary microseconds (4096 seconds)
;   There can be a time skip after the first <~1000 events (the first packet?)
;     This may result from hitting "Save to File" before "AcquireRT"

t=d.time/(2d^20)
;n = n_elements(t)

skip = where((t-shift(t,1))[1:*] gt 1,nskip)+1
;Skip correction is restricted to a single, >1 second skip within the first 1000 events
if nskip eq 1 and skip le 1000 then begin
  d=d[skip:*]
  t=t[skip:*]-t[skip[0]]
endif

k = where((t-shift(t,1))[1:*] lt -2d^11,nk)+1
if nk gt 0 then for i=0,nk-1 do t[k[i]:*] += 2d^12

; Old method of calculating time
;timearr = double(d[where(d.channel ne 0)].time)
;time = timearr[n_elements(timearr)-1]-timearr[0]+total((timearr-shift(timearr,1))[1:*] lt -2d^31)*2d^32
;time /= 2d^20

time = max(t)-min(t)

front = where(d.a2d_index eq 8, nfront)
rear = where(d.a2d_index eq 17, nrear)

tbin = fcheck(slice, time+2d^(-20))

nslice = ceil(time/tbin)

if nfront gt 0 then begin
  lc_front = histogram(t[front],min=0,bin=1)
  spec_front = hist_2d(d[front].channel,t[front],min1=-2,max1=8191,min2=0,max2=time,bin2=tbin)
endif else begin
  lc_front = lonarr(ceil(time))
  spec_front = reform(lonarr(8194,nslice))
endelse
uld_front = spec_front[1,*] ; -1
reset_front = spec_front[0,*] ; -2
if reset_front gt 0 then begin
  temp_use = where(d[front].channel eq -2)
  temp_offset = ((d[front])[temp_use].time-shift((d[front])[temp_use].time,1))[1:*]
  ;print,(histogram(temp_offset,bin=10))[0:10]
  reset_front = total(temp_offset ge 40) ; ad-hoc value to call a reset a duplicate
endif
spec_front = spec_front[2:*,*]

if nrear gt 0 then begin
  lc_rear = histogram(t[rear],min=0,bin=1)
  spec_rear = hist_2d(d[rear].channel,t[rear],min1=-2,max1=8191,min2=0,max2=time,bin2=tbin)
endif else begin
  lc_rear = lonarr(ceil(time))
  spec_rear = reform(lonarr(8194,nslice))
endelse
uld_rear = spec_rear[1,*] ; -1
reset_rear = spec_rear[0,*] ; -2
if reset_front gt 0 then begin
  temp_use = where(d[rear].channel eq -2)
  temp_offset = ((d[rear])[temp_use].time-shift((d[rear])[temp_use].time,1))[1:*]
  ;print,(histogram(temp_offset,bin=10))[0:10]
  reset_rear = total(temp_offset ge 40) ; ad-hoc value to call a reset a duplicate
endif
spec_rear = spec_rear[2:*,*]

if keyword_set(plot) then begin
  temp = !p.multi

  !p.multi = [0,1,nslice+1]

  maxcounts = max(spec_front) > max(spec_rear)

  for i = 0,nslice-1 do begin
    plot,spec_front[*,i],yr=[1,maxcounts],/yl,psym=10,$
      xtitle='Channel',ytitle='Counts',title='Count spectrum'
    oplot,spec_rear[*,i],psym=10,color=6
  endfor

  plot,lc_front,yr=[1,max([lc_rear,lc_front])],psym=10,$
    xtitle='Seconds',ytitle='Total counts',title='Count lightcurve'
  oplot,lc_rear,psym=10,color=6
  for i=0,nslice-1 do oplot,(i ne nslice-1 ? tbin*(i+1) : time)*[1,1],!y.crange,color=5,linestyle=1

  !p.multi = temp
endif

obj_destroy,o

if keyword_set(save) then begin
  stub=strsplit(filename,'.ssr',/extract)
  save,filename,spec_front,spec_rear,lc_front,lc_rear,slice,time,filename=stub+'.sav'
endif

tab = string(9b)

print,"*** Front segment ***"
for i = 0,nslice-1 do begin
  duration = float(i ne nslice-1 ? tbin : time-(nslice-1)*tbin)
  print,"Duration: "+tab,duration," s"
  print,"Count rate: "+tab,total(spec_front[*,i])/duration," counts/s"
  print,"Reset period: "+tab,duration/reset_front[i]," s/reset"
  print,"Reset rate: "+tab,reset_front[i]/duration," resets/s"
  print,"ULD rate: "+tab,uld_front[i]/duration," ULDs/s"
  if keyword_set(source) then autoroi,spec_front[*,i]/duration,source,gain=gain_front
endfor

print,"*** Rear segment ***"
for i = 0,nslice-1 do begin
  duration = float(i ne nslice-1 ? tbin : time-(nslice-1)*tbin)
  print,"Duration: "+tab,duration," s"
  print,"Count rate: "+tab,total(spec_rear[*,i])/duration," counts/s"
  print,"Reset period: "+tab,duration/reset_rear[i]," s/reset"
  print,"Reset rate: "+tab,reset_rear[i]/duration," resets/s"
  print,"ULD rate: "+tab,uld_rear[i]/duration," ULDs/s"
  if keyword_set(source) then autoroi,spec_rear[*,i]/duration,source,/rear,gain=gain_rear
endfor

end
