pro check_SED

restore, 'ascii_templates/margestats_template.sav'
file_low = findfile('SED_fitting/COSMOS_SEDfits/o2sed_40kmet_lowprior/*.margestats')
file_hi = findfile('SED_fitting/COSMOS_SEDfits/o2sed_40kmet/*.margestats')

age_low = dblarr(n_elements(file_low))
age_hi = dblarr(n_elements(file_hi))                                                                        

for i = 0, n_elements(file_low)-1 do begin
   marge_low = read_ascii(file_low[i], template = margestats_template)
   age_low[i] = marge_low.field2[0]
endfor

for i = 0, n_elements(file_hi)-1 do begin
   marge_hi = read_ascii(file_hi[i], template = margestats_template)
   age_hi[i] = marge_hi.field2[0]
endfor
                                                      
age_low = age_low*alog10(2.71828182846)      ; log(age)            
age_hi = age_hi*alog10(2.71828182846)

ind_low = fix(strmid(file_low, 52, 3))  ;Change '52' depending on lenght of file names!!!!!                            
ind_hi = fix(strmid(file_hi, 43,3))

b = where(age_hi lt 6)      ; These are the ones of the initial SED fit that are hitting the prior
c = ind_hi[b]               ; These are the indices of the ages hitting the prior
                                ; ind_low are the SEDs that managed to
                                ; come into existence after changing
                                ; the prior

stop

end
