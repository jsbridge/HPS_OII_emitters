;+
;
; Looks at CAS, Gini, M20, Sersic indices of [OII] emitters in COSMOS
; data
; Note: Zurich cataloge doesn't have ACS tile # column as shown
; in README on COSMOS site, so all field numbers are shifted up one
; after field03
; 
;-
pro OII_morph

;restore, 'ascii_templates/zamoj_template.sav'
;morph_data = read_ascii('data/cosmos_morph_zamoj_HPSfield.tbl', template = zamoj_template)

restore, 'ascii_templates/zurich_template.sav'
morph_data = read_ascii('data/cosmos_morph_zurich_HPSfield.tbl', template = zurich_template)


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

;plot, morph_data.field02[cos_min], morph_data.field03[cos_min], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
;oplot, HPSRA[140:323], HPSdec[140:323], psym=4, color=240

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
stop
;plot, morph_data.field02[cos_min], morph_data.field03[cos_min], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
;oplot, cos_OIIRA, cos_OIIdec, psym=4, color=240

cosmorphRA = morph_data.field02[cos_min]
cosmorphdec = morph_data.field03[cos_min]

;plot, cosmorphRA[cos_OII - 140], cosmorphdec[cos_OII - 140], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
;oplot, OIIRA, OIIdec, psym=4, color=240

cos_sersic = (morph_data.field56[cos_min])[cos_OII - 140]
cos_index = where(cos_sersic eq 'null')
remove, cos_index, cos_sersic

cos_conc = (morph_data.field31[cos_min])[cos_OII - 140]
cos_asym = (morph_data.field32[cos_min])[cos_OII - 140]
cos_M20 = (morph_data.field30[cos_min])[cos_OII - 140]
cos_Gini = (morph_data.field29[cos_min])[cos_OII - 140]

openpps, 'morphology/cosmos_sersic'
plothist, double(cos_sersic), bin = 0.15, title = 'Sersic Indices of [OII] Emitters in COSMOS field of HPS', xtitle = 'Sersic Index', ytitle = 'Number'
closepps

openpps, 'morphology/cosmos_asym'
plothist, cos_asym, bin = 0.01, title = 'Asymmetry of [OII] Emitters in COSMOS field of HPS', xtitle = 'A', ytitle = 'Number'
closepps

openpps, 'morphology/cosmos_conc'
plothist, cos_conc, bin = 0.1, title = 'Concentration of [OII] emitters in COSMOS field of HPS', xtitle = 'C', ytitle = 'Number'
closepps

openpps, 'morphology/cosmos_M20_vs_Gini'
plot, cos_M20, cos_Gini, yrange=[0.3, 0.75], ystyle=1, xrange=[-.4, -3], xstyle=1, psym = 6, title = 'M20 vs. Gini Coefficient for [OII] emitters in COSMOS field of HPS', xtitle = textoidl('M_{20}'), ytitle = 'G'
oplot, [-0, -3], [0.375, 0.708], linestyle=2 ; line for normal galaxies according to Lotz et al.
closepps

openpps, 'morphology/cosmos_M20_vs_conc'
plot, cos_M20, cos_conc, psym=6, yrange = [1, 5], ystyle=1, xrange =[-.4, -3], xstyle=1, title = 'M20 vs. Concentration for [OII] emitters in COSMOS field of HPS', xtitle = textoidl('M_{20}'), ytitle = 'C'
closepps

openpps, 'morphology/cosmos_asym_vs_conc'
plot, cos_asym, cos_conc, psym=6, yrange = [1, 5], ystyle=1, xrange =[0, 0.6], xstyle=1, title = 'Asymmetry vs. Concentration for [OII] emitters in COSMOS field of HPS', xtitle = 'A', ytitle = 'C'
closepps

openpps, 'morphology/fullcosmos_M20_vs_Gini'
plot, morph_data.field30, morph_data.field29, psym = 6, yrange=[0.3, 0.75], ystyle=1, xrange=[-.4, -3], xstyle=1, title = 'M20 vs. Gini Coefficient for [OII] emitters in COSMOS field of HPS (blue)', xtitle = textoidl('M_{20}'), ytitle = 'G', charsize=1 ; ALL COSMOS points
oplot, cos_M20, cos_Gini, psym = 4, color=60
closepps

openpps, 'morphology/fullcosmos_M20_vs_conc'
plot, morph_data.field30, morph_data.field31, psym=6, yrange = [1, 5], ystyle=1, xrange =[-.4, -3], xstyle=1, title = 'M20 vs. Concentration for [OII] emitters in COSMOS field of HPS (red)', xtitle = textoidl('M_{20}'), ytitle = 'C', charsize=1
oplot, cos_M20, cos_conc, psym=4, color=240
closepps

openpps, 'morphology/fullcosmos_asym_vs_conc'
plot, morph_data.field32, morph_data.field31, psym=6, yrange = [1, 5], ystyle=1, xrange =[0, 0.6], xstyle=1, title = 'Asymmetry vs. Concentration for [OII] emitters in COSMOS field of HPS (green)', xtitle = 'A', ytitle = 'C', charsize=1
oplot, cos_asym, cos_conc, psym=4, color=150
closepps

;remove, cos_index, cos_asym, cos_conc, cos_M20, cos_Gini
;cos_sm_sers_M20 = dblarr(100)
;cos_sm_sers_Gini = dblarr(100)
;j = 0
;for i = 0, n_elements(cos_sersic)-1 do begin
;   if (double(cos_sersic))[i] le 1 then begin
;      cos_sm_sers_M20[j] = cos_M20[i]
;      cos_sm_sers_Gini[j] = cos_Gini[i]
;      j += 1
;   endif
;endfor

;openpps, 'morphology/cos_sm_sersic_M20_vs_Gini'
;plot, cos_M20, cos_Gini, yrange=[0.3, 0.75], ystyle=1, xrange=[-.4, -3], xstyle=1, psym = 6, title = 'M20 vs. Gini Coefficient for [OII] emitters in COSMOS field of HPS, Sersic indices < 0.5 in green', charsize=.9, xtitle = textoidl('M_{20}'), ytitle = 'G'
;oplot, cos_sm_sers_M20, cos_sm_sers_Gini, psym = 4, color=150, thick=7
;closepps

cos_halfrad = (morph_data.field41[cos_min])[cos_OII - 140]
cos_index = where(cos_halfrad eq 'null')
remove, cos_index, cos_halfrad, cos_Gini, cos_M20, cos_conc, cos_asym

openpps, 'morphology/cosmos_halfrad_vs_gini'
plot, cos_halfrad, cos_Gini, psym = 6, title = 'Gini coefficient vs. Half-light Radius for [OII] Emitters in COSMOS field of HPS', xtitle = 'Half-light radius (arcsec)', ytitle = 'G', charsize = 1
closepps

openpps, 'morphology/cosmos_halfrad_vs_M20'
plot, cos_halfrad, cos_M20, psym = 6,  title = 'M20 vs. Half-light Radius for [OII] Emitters in COSMOS field of HPS', xtitle = 'Half-light radius (arcsec)', ytitle = textoidl('M_{20}')
closepps

openpps, 'morphology/cosmos_halfrad_vs_conc'
plot, cos_halfrad, cos_conc, psym = 6, title = 'Concentration vs. Half-light Radius for [OII] Emitters in COSMOS field of HPS', xtitle = 'Half-light radius (arcsec)', ytitle = 'C', charsize = 1
closepps

openpps, 'morphology/cosmos_halfrad_vs_asym'
plot, cos_halfrad, cos_asym, psym = 6, title = 'Asymmetry vs. Half-light Radius for [OII] Emitters in COSMOS field of HPS', xtitle = 'Half-light radius (arcsec)', ytitle = 'A', charsize = 1
closepps

openpps, 'morphology/cosmos_halfrad'
plothist, cos_halfrad, bin = .15, title = 'PSF-convolved Half-light Radii of [OII] Emitters in COSMOS field of HPS', xtitle = 'Half-light radius (arcsec)', ytitle = 'Number'
closepps

; And onto GOODS-N!

goods_CAS = mrdfits('data/GOODS_N_morphs.fits',1) ; Not great matching :(
goods_sers = mrdfits('data/ACS-GC_published_catalogs/goods_v_i_public_catalog_V1.0.fits',1)
;restore, 'ascii_templates/goods_morph_template.sav'
;goods_nic = read_ascii('data/goodsn_nicmos_morph.txt', template = goods_morph_template) ; no match up whatsoever, oh well

gval1_sers = dblarr(144, n_elements(goods_sers.RA))
gval2_sers = dblarr(144, n_elements(goods_sers.RA))
gval1_CAS = dblarr(144, n_elements(goods_CAS.RA))
gval2_CAS = dblarr(144, n_elements(goods_CAS.RA))
;gval1_nic = dblarr(144, n_elements(goods_nic. field01))
;gval2_nic = dblarr(144, n_elements(goods_nic.field01))
j=0
for i = 324, 467 do begin
   gval1_sers[j, *] = goods_sers.RA - HPSRA[i]
   gval2_sers[j, *]= goods_sers.dec - HPSdec[i]
   gval1_CAS[j, *] = goods_CAS.RA - HPSRA[i]
   gval2_CAS[j, *]= goods_CAS.dec - HPSdec[i]
;   gval1_nic[j, *] = goods_nic.field02 - HPSRA[i]
;   gval2_nic[j, *]= goods_nic.field03 - HPSdec[i]
   j+= 1
endfor

gdist_sers = dblarr(144, n_elements(goods_sers.RA))
gdist_CAS = dblarr(144, n_elements(goods_CAS.RA))
;gdist_nic = dblarr(144, n_elements(goods_nic.field01))
for i = 0,143 do begin
   gdist_sers[i, *] = sqrt(cos(gval2_sers[i,*] * !pi/180d)^2 * gval1_sers[i, *]^2 + gval2_sers[i, *]^2)
   gdist_CAS[i, *] = sqrt(cos(gval2_CAS[i,*] * !pi/180d)^2 * gval1_CAS[i, *]^2 + gval2_CAS[i, *]^2)
;   gdist_nic[i, *] = sqrt(cos(gval2_nic[i,*] * !pi/180d)^2 * gval1_nic[i, *]^2 + gval2_nic[i, *]^2)
endfor

goods_min_sers = intarr(144)
goods_min_CAS = intarr(144)
;goods_min_nic = intarr(144)
for i = 0, 143 do begin
   goods_min_sers[i] = where(gdist_sers[i,*] eq min(gdist_sers[i, *]))   ; returns indices of cos_data array of the likely objects
   goods_min_CAS[i] = where(gdist_CAS[i,*] eq min(gdist_CAS[i, *]))
;   goods_min_nic[i] = where(gdist_nic[i,*] eq min(gdist_nic[i, *]))
endfor

goodsCASRA = goods_CAS.RA
goodsCASdec = goods_CAS.dec
goods_sersRA = goods_sers.RA
goods_sersdec = goods_sers.dec
goods_sersic = goods_sers.n_galfit_hi

goods_OII = O_index[where(O_index ge 324 and O_index le 467)]
goods_sersic = (goods_sersic[goods_min_sers])[goods_OII - 324]
goods_RA = HPSRA[324:467]
goods_dec = HPSdec[324:467]
goods_OIIRA = goods_RA[goods_OII - 324]
goods_OIIdec = goods_dec[goods_OII - 324]

;window, 0
;plot, goodsCASRA[goods_min_CAS], goodsCASdec[goods_min_CAS], psym=6, title='GOODS-N (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[62.17, 62.26],ystyle=1, xrange=[189.05, 189.3], xstyle=1
;oplot, HPSRA[324:467], HPSdec[324:467], psym=4, color=240

;window, 1
;plot, (goods_sersRA[goods_min_sers])[goods_OII - 324], (goods_sersdec[goods_min_sers])[goods_OII - 324], psym=6, title='GOODS-N (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[62.17, 62.26], ystyle=1, xrange=[189.05, 189.3], xstyle=1
;oplot, goods_OIIRA, goods_OIIdec, psym=4, color=240

;window, 2
;plot, goods_nic.field02[goods_min_nic], goods_nic.field03[goods_min_nic], psym=6, title='GOODS-N (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[62.17, 62.26], ystyle=1, xrange=[189.05, 189.3], xstyle=1
;oplot, HPSRA[324:467], HPSdec[324:467], psym=4, color=240

loadct, 39
openpps, 'morphology/goods_sersic'
plothist, goods_sersic, bin = 0.15, title = 'Sersic Indices of [OII] Emitters in GOODS-N field of HPS', xtitle = 'Sersic Index', ytitle = 'Number'
closepps
openpps, 'morphology/goods_cos_sersic'
plothist, goods_sersic, bin = 0.15, title = 'Sersic Indices of [OII] Emitters in GOODS-N (black) and COSMOS (red) field of HPS', xtitle = 'Sersic Index', ytitle = 'Number', yrange = [0, 25], ystyle=1, xrange = [0, 5.6], xstyle=1, charsize=1
plothist, double(cos_sersic), bin = 0.15, /overplot, color = 240
closepps
stop
end
