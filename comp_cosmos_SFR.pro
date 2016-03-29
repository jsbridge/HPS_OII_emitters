pro comp_cosmos_SFR

cosmos_OII = cosmos_OII()

restore, 'ascii_templates/margestats_template.sav'

mass = dblarr(n_elements(cosmos_OII.id))
SFR = dblarr(n_elements(cosmos_OII.id))
age = dblarr(n_elements(cosmos_OII.id))

; This program calls the SED fits that had 40k runs, and fit for
; metallicity, age, stellar mass, galaxy mass, and E(B-V).
; It assumes a constant star formation.

for i = 0, n_elements(cosmos_OII.id)-1 do begin
   marge = read_ascii('SED_fitting/COSMOS_SEDfits/o2sed_40kmet/hps'+strcompress(string(cosmos_OII.id[i]), /remove_all)+'.margestats', template = margestats_template)
   mass[i] = marge.field2[4]
   SFR[i] = marge.field2[4]/marge.field2[0]
   age[i] = marge.field2[0]
endfor

mass = mass*alog10(2.71828182846)   ; log(mass [M_sol])
SFR = SFR*alog10(2.71828182846)       ; log(SFR [M_sol/yr])
SSFR = alog10(10^(SFR)/10^(mass))     ; log(SSFR [yr^-1])
age = age*alog10(2.71828182846)      ; log(age)

readcol, 'data/robin_sfr.dat', robin_id, robin_RA, robin_dec, robin_UV, robin_OII, format = 'I, D, D, D, D'

arr = intarr(n_elements(cosmos_OII.id))
for i = 0, n_elements(cosmos_OII.id) - 1 do begin
   arr[i] = findel(robin_id, cosmos_OII.id[i])
endfor

loadct, 39
openpps, 'SFR_plots/comp_cosmos_SFR_OII'
plot, mass, robin_OII[arr], psym = 6,  xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-1.5, 2.1]
oplot, mass, SFR, psym = 6, color = 200
legend, ['SED fitting SFR', '[OII] Line SFR'], color = [200, 0], psym = [6, 6]
closepps
openpps, 'SFR_plots/comp_cosmos_SFR_UV'
plot, mass, robin_UV[arr], psym = 6,  xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-.05, .8]
oplot, mass, SFR, psym = 6, color = 100
legend, ['SED fitting SFR', 'UV SFR'], color = [100, 0], psym = [6, 6]
closepps
openpps, 'SFR_plots/cosmos_scatt_SFR_OII'
plot, SFR, robin_OII[arr], psym =6, xtitle = textoidl('log(SED [OII] SFR [M_{sun} yr^{-1}])'), ytitle = textoidl('log([OII] Line SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-2, 2.2]
closepps
openpps, 'SFR_plots/cosmos_scatt_SFR_UV'
plot, SFR, robin_UV[arr], psym =6, xtitle = textoidl('log(SED [OII] SFR [M_{sun} yr^{-1}])'), ytitle = textoidl('log(UV SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-.02, 0.5]
closepps
openpps, 'SFR_plots/SEDtime_vs_OIIlineSFR'
plot, robin_OII[arr], age, psym =6, xtitle = textoidl('log([OII] line SFR [M_{sun} yr^{-1}])'), ytitle = textoidl('log(SED-derived Age [yrs])'), ystyle=1, yrange=[5, 10]
closepps
openpps, 'SFR_plots/SEDtime_vs_OIIlineSSFR'
plot, robin_OII[arr]/mass, age, psym =6, xtitle = textoidl('log([OII] line SSFR [yr^{-1}])'), ytitle = textoidl('log(SED-derived Age [yrs])') , ystyle=1, yrange=[5, 10]
closepps

stop
end
