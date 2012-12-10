; NAME: getgain
;
; PURPOSE: Produces a energy gain calibration using a linear fit to several
;	points of channel and corresponding energy (e.g., from "roi").  It is
;	significantly easier to use "autoroi" if you can.
;
; INPUTS:
;   channels	list of channels
;   energies	list of corresponding energies
;   plot	keyword - plots residuals of linear fit
;
; OUTPUTS:
;   (return)	1D energy gain calibration
;   coeff	[Y-intercept, slope] from linear fit
;
; EXAMPLES:
;   gf = getgain([183.4,262.3,818.3,893.4,1044.4,1123.6],[53.2,81.0,276.4,302.9,356.0,383.8])
;
; HISTORY:
;   2012-12-08, AYS: release

function getgain,channels,energies,plot=plot,coeff=coeff

coeff = linfit(channels,energies,chisq=chisq)

print,coeff,chisq,format='("E=",F7.2,"+",F7.5,"*C, chisq=",F7.4)'

x = findgen(8192)
y = coeff[0]+coeff[1]*findgen(8192)

;plot,x,y
;oplot,channels,energies,psym=5,color=6

n = n_elements(energies)

if keyword_set(plot) then plot,indgen(n)+1,coeff[0]+channels*coeff[1]-energies,psym=4,xr=[0,n+1],/xs,xticks=n+1,xtickname=[0,energies,0],xtitle='Energy (keV)',ytitle='Error (keV)'

return,y

end
