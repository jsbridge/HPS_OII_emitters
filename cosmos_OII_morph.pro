pro cosmos_OII_morph

cat = mrdfits('data/ACS-GC_published_catalogs/cosmos_i_public_catalog_V1.0.fits.gz',1)

restore, 'ascii_templates/cassata_template.sav'
morph_data = read_ascii('data/cosmos_morph_cassata_1.1.tbl', template = cassata_template)

; First thing, pull out all galaxies from cat that have photo z of
; <0.57

cosRA = cat.RA
cosdec = cat.dec
cosmos_z = cat.photoz

ind = where((cosmos_z lt 0.0) or (cosmos_z gt 0.57))
remove, ind, cosRA, cosdec, cosmos_z

; Now that I have all the galaxies with z < 0.57, I need to correspond
; those to galaxies in the morphology catalog
ind = dblarr(n_elements(cosRA))
for i = 0, n_elements(cosRA) - 1 do begin
   ind[i] = findel(cosRA[i], morph_data.field02)
endfor

allcos_propz_M20 = morph_data.field10[ind]
allcos_propz_Gini = morph_data.field09[ind]


; The following code pulls out the morphologies of the OII emitters

readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec

val1 = dblarr(184, n_elements(morph_data.field02))
val2 = dblarr(184, n_elements(morph_data.field03))
j=0
for i = 140, 323 do begin
   val1[j, *] = morph_data.field02 - HPSRA[i]
   val2[j, *] = morph_data.field03 - HPSdec[i]
   j+= 1
endfor

dist = dblarr(184, n_elements(morph_data.field02))
for i = 0,183 do begin
   dist[i, *] = sqrt(cos(val2[i,*] * !pi/180d)^2 * val1[i, *]^2 + val2[i, *]^2)
endfor

cos_min = lonarr(184)
min = dblarr(184)
for i = 0, 183 do begin
   cos_min[i] = where(dist[i,*] eq min(dist[i, *]))   ; returns indices of cos_data array of the likely objects
   min[i] = min(dist[i, *])
endfor

restore, 'ascii_templates/HPS4_template.sav'
HPS_data = read_ascii('data/HPStable4.txt', template=template)

oxy = HPS_data.field11
remove, ind, oxy

O_index = where(oxy eq '[OII]')
cos_OII = O_index[where(O_index ge 140 and O_index le 323)]
cos_RA = HPSRA[140:323]
cos_dec = HPSdec[140:323]
cos_OIIRA = cos_RA[cos_OII - 140]
cos_OIIdec = cos_dec[cos_OII - 140]

cosmorphRA = morph_data.field02[cos_min]
cosmorphdec = morph_data.field03[cos_min]

;plot, cosmorphRA[cos_OII - 140], cosmorphdec[cos_OII - 140], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
;oplot, cos_OIIRA, cos_OIIdec, psym=4, color=240

cos_OII_M20 = (morph_data.field10[cos_min])[cos_OII - 140]
cos_OII_Gini = (morph_data.field09[cos_min])[cos_OII - 140]
cos_halflight = (morph_data.field06[cos_min])[cos_OII - 140]
stop
; Plotting all COSMOS galaxies > 0.57 as contours and then overplotting [OII]
; galaxies 
; First get rid of stupid Gini or M20 values (like outside a
; reasonable range) else the contours freak out... I think.
new_Gini = allcos_propz_Gini
new_M20 = allcos_propz_M20
ind = where(new_Gini gt 0.85 or new_Gini lt 0.05)
remove, ind, new_Gini, new_M20 
ind = where(new_M20 gt 0 or new_M20 lt -3)
remove, ind, new_gini, new_M20

loadct, 39
openpps, 'SFR_plots/M20vsG_contours'
dens2d, new_M20, new_Gini, g, xtitle = textoidl('M_{20}'), ytitle = 'Gini' , psym = -1, /contur, nbin = 60, charsize = 1.2;, yrange=[0.05, 0.85], xrange=[0, -3]
oplot, cos_OII_M20, cos_OII_Gini, psym = 6, color = 245
oplot, [0, -3], [0.33, 0.75], linestyle=2
closepps

; The following chunk of code gets the masses of the OII emitters
cosmos_OII = cosmos_OII()

restore, 'ascii_templates/likestats_template_tau.sav'   ; will need to change template if not doing exponential SFH
file = findfile('SED_fitting/COSMOS_SEDfits/o2sed_40ktau/*.likestats')  ; modify the folder name for the different fits

mass = dblarr(n_elements(file))

; This program calls the SED fits.  I've just been modifying it
; to work for various versions and parameters

for i = 0, n_elements(file)-1 do begin
   data = read_ascii(file[i], template = likestats_template_tau)   ;  Check you've got the right template
   mass[i] = data.field2[4]    ; Check the subscripting numbers here - depending on what params have been fit, they may be in different rows or whatnot
endfor

mass = mass*alog10(2.71828182846)   ; log(mass [M_sol])

openpps, 'SFR_plots/cosmos_M20vsmass'
plot, mass, cos_OII_M20, psym = 6, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('M_{20}')
oplot, mass, cos_OII_M20, psym = 6, color = 60
closepps

openpps, 'SFR_plots/cosmos_Ginivsmass'
plot, mass, cos_OII_Gini, psym = 6,  xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = 'Gini'
oplot, mass, cos_OII_Gini, psym = 6, color = 60
closepps

print, 'M20 vs. Mass: r=', correlate(mass, cos_OII_M20)
print, 'Gini vs. Mass: r=', correlate(mass, cos_OII_Gini)

stop
end
