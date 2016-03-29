pro cosmos_mass_hist

cosmos_OII = cosmos_OII()

restore, 'ascii_templates/margestats_template.sav'

mass = dblarr(n_elements(cosmos_OII.id))

for i = 0, n_elements(cosmos_OII.id)-1 do begin
   marge = read_ascii('SED_fitting/COSMOS_SEDfits/o2sed_40kmet/hps'+strcompress(string(cosmos_OII.id[i]), /remove_all)+'.margestats', template = margestats_template)
   mass[i] = marge.field2[4]
endfor

mass = mass*alog10(2.71828182846)   ; log(mass [M_sol])

print, n_elements(mass)

openpps, 'cosmos_mass_hist'
plothist, mass, /autobin, xtitle = textoidl('log(mass [M_{sol}])')
closepps

end
