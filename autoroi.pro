; NAME: autoroi
;
; PURPOSE: Uses "good" channels to auto-calibrate a spectrum based on a
;	specified radioactive source
;
; INPUTS:
;   spectrum	1D spectrum
;   source	integer, positive for front segment, negative for rear segment
;		1: Ba-133
;		2: Co-60/Eu-152/Eu-154
;   rear	keyword - if set, forces negative sign for source
;   plot	keyword - passthrough to "getgain", plots gain calibration
;
; OUTPUTS:
;   result	array with fitted line energies, FWHMs, and strengths
;   gain	energy gain calibration
;
; EXAMPLES:
;   autoroi,sf,1
;   autoroi,sr,-2
;
; HISTORY:
;   2012-12-08, AYS: release

pro autoroi,spectrum,source,result=result,gain=gain,plot=plot,rear=rear

source2 = keyword_set(rear) ? -abs(source) : source

case source2 of
  1: begin
      band = [[175,190],[250,275],[810,835],[880,915],[1020,1065],[1110,1145]]
      energies = [53.1625,80.9971,276.3997,302.8510,356.0134,383.8480]
    end
  -1: begin
      band = [[215,240],[715,750],[780,820],[915,955],[990,1025]]
      energies = [80.9971,276.3997,302.8510,356.0134,383.8480]
    end
  2: begin
      band = [[2070,2105],[3340,3390],[3635,3680],[3800,3855],[4010,4065]]
      energies = [723.305,1173.237,1274.436,1332.501,1408.006]
    end
  -2: begin
      band = [[1860,1905],[3010,3075],[3275,3330],[3415,3480],[3615,3675]]
      energies = [723.305,1173.237,1274.436,1332.501,1408.006]
    end
endcase

roi,spectrum,band,result=result

gain = getgain(reform(result[0,*]),energies,plot=plot,coeff=coeff)

result *= coeff[1]
result[0,*] += coeff[0]

print,result

end
