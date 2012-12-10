
pro fu_gauss_lin, x,a,f,pder
ex = exp(-(x-a(2))^2./(2.*a(1)^2))
exm = a(0)/sqrt(2.*!pi)/a(1)*ex
f = exm  + a(4)*x + a(3)

if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
pder(*,0) = exm/a(0)
pder(*,1) = -exm/a(1) + exm*(x-a(2))^2./a(1)^3.
pder(*,2) = exm*(x-a(2))/a(1)^2.
pder(*,3) = 1.
pder(*,4) = x
end
