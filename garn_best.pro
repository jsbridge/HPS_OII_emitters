function garn_best, mass

k_Ha = 2.659*(-1.857+1.04/0.6563)+4.05

E_BV = dblarr(n_elements(mass))

for i = 0, n_elements(mass)-1 do begin
    x = alog10((10^(mass[i]))/(1d10))
    E_BV[i] = (0.91 + 0.77*x + 0.11*x^2 - 0.09*x^3)/k_Ha
endfor

return, E_BV

end
