pro fitgauss2,x=x,y,w=w,fit_fwhm,fit_center,fit_area,ffe,fce,fae,quiet=quiet,$
   startwidth=startwidth, startcent=startcent,yfit=yfit,chisq=chisq,estimates

x=fcheck(x,findgen(n_elements(y)))
w=fcheck(w,fltarr(n_elements(y))+1.)
y1=y[0]
y2=y[n_elements(y)-1]
x1=x[0]
x2=x[n_elements(x)-1]
baseave=(y1+y2)/2.
baseslope=(y2-y1)/(x2-x1)
center = (where(y EQ max(y)))[0]
baseline = baseave + baseslope*(x-(x[center])[0])
height=max(y-baseline)
peak = where(y-baseline GT height/2.)

width = (x[max(peak)]-x[min(peak)]) / (2.*sqrt(2.*alog(2.)))
if width EQ 0 then width =( x[1]-x[0] ) / (2.*sqrt(2.*alog(2.)))

estimates=[height*width*sqrt(2.*!pi),width,x[center],$
           baseave-baseslope*center,baseslope]
if (keyword_set(startwidth)) then estimates[1]=startwidth
if (keyword_set(startcent)) then estimates[2]=startcent

if NOT keyword_set(quiet) then print,estimates
for j=0, 5 do begin
     yfit=curvefit(x,y,w,estimates,chisq=chisq,sigma,/double,function_name='fu_gauss_lin')
end

if NOT keyword_set(quiet) then begin
  plot,x,y
  oploterr,x,y,sqrt(1./w)
  oplot,x,yfit,psym=0
;  print,estimates
;  print,sigma
endif
fit_fwhm = estimates[1]*2.*sqrt(2.*alog(2.))
fit_center = estimates[2]
fit_area = estimates[0]
ffe = sigma[1]*2.*sqrt(2.*alog(2.))
fce = sigma[2]
fae = sigma[0]
end


