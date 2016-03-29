;+
;
; Compares spec z's of COSMOS and GOODS-N to HPS redshifts
; (except only some specs z's of COSMOS are available, so the
; rest are photo z's)
;
;-
pro comp_zspec

restore, 'ascii_templates/HPS4_template.sav'
data = read_ascii('data/HPStable4.txt', template=template)
z = data.field12

; This section is to acquire what spectroscopic redshifts are
; available for COSMOS
readcol, 'data/HSTCOSMOS_ZCOSBRIGHTSPEC10k.txt', ID, zspec, specRA, specdec, format='L, D, X, X, D, D', comment='#'
readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec, z

val1_spec = dblarr(184, n_elements(specRA))
val2_spec = dblarr(184, n_elements(specdec))
j=0
for i = 140, 323 do begin
   val1_spec[j, *] = specRA - HPSRA[i]
   val2_spec[j, *]= specdec - HPSdec[i]
   j+= 1
endfor

dist_spec = dblarr(184, n_elements(specRA))
for i = 0,183 do begin
   dist_spec[i, *] = sqrt(cos(val2_spec[i,*] * !pi/180d)^2 * val1_spec[i, *]^2 + val2_spec[i, *]^2)
endfor

cos_min_spec = intarr(184)
min_spec = dblarr(184)
for i = 0, 183 do begin
   cos_min_spec[i] = where(dist_spec[i,*] eq min(dist_spec[i, *]))   ; returns indices of cos_data array of the likely objects
   min_spec[i] = min(dist_spec[i, *])
endfor

;newRA = specRA[cos_min_spec]
;newdec = specdec[cos_min_spec]

get_zspec = where(min_spec lt 5d-4, complement = get_zphot)    ; This is to only grab the spec z's that exist for the HPS sources, get_zphot are the indices must use photo z's instead

;plot, newRA[get_zspec], newdec[get_zspec], psym = 6
;oplot, HPSRA[140:323], HPSdec[140:323], psym=4, color=240

loadct, 39
openpps, 'z_comp/corresp_pts_COSMOS_zspec'   
plot, specRA[cos_min_spec], specdec[cos_min_spec], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
oplot, HPSRA[140:323], HPSdec[140:323], psym=4, color=240
closepps

; This section is to get the photometric redshifts for the
; things we don't have spectroscopic redshifts for from COSMOS
restore, 'ascii_templates/cosmos_template.sav'
cos_data = read_ascii('data/cosmos_zphot_HPSfield.tbl', template = cosmos_template)

val1_phot = dblarr(184, n_elements(cos_data.field03))
val2_phot = dblarr(184, n_elements(cos_data.field04))
j=0
for i = 140, 323 do begin
   val1_phot[j, *] = cos_data.field03 - HPSRA[i]
   val2_phot[j, *]= cos_data.field04 - HPSdec[i]
   j+= 1
endfor

dist_phot = dblarr(184, n_elements(cos_data.field03))
for i = 0,183 do begin
   dist_phot[i, *] = sqrt(cos(val2_phot[i,*] * !pi/180d)^2 * val1_phot[i, *]^2 + val2_phot[i, *]^2)
endfor

cos_min_phot = intarr(184)
for i = 0, 183 do begin
   cos_min_phot[i] = where(dist_phot[i,*] eq min(dist_phot[i, *]))
endfor

; want subset of cos_min_spec... subscript by get_zspec, then I will
; have indices of initial set that qualify, so then subscript zspec by
; that subscription
yay = cos_min_spec[get_zspec]
yayz = zspec[yay]
cosmosz = dblarr(184)
cosmosz[get_zspec] = yayz
boo = cos_min_phot[get_zphot]
booz = cos_data.field05[boo]
cosmosz[get_zphot] = booz

; Okay, but now, how do I mine a complete different data set with
; different indices for the photo z's, and how do I keep the
; redshifts straight with which z belongs to which object?
; Well, I've got the indices of the HPS objects that are covered with get_zspec

cosmos_z = cosmosz
new1_z =  z[140:323]
yay = where(cosmos_z gt 0 and cosmos_z lt 5, complement = bad)
remove, bad, cosmos_z, new1_z 

openpps, 'z_comp/HPS_COSMOS_zspecphot'
plot, cosmosz, z[140:323], psym=6, xrange=[0, 4], xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS spec z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(cosmos_z, new1_z)))
closepps

; This section is to get GOODS-N spectroscopic redshifts
restore, 'ascii_templates/goods_temp.sav'
goods_data = read_ascii('data/goodsn_multiwave_v2.0.cat', template = goods_template)

gval1 = dblarr(144, n_elements(goods_data.field01))
gval2 = dblarr(144, n_elements(goods_data.field01))
j=0
for i = 324, 467 do begin
   gval1[j, *] = goods_data.field03 - HPSRA[i]
   gval2[j, *]= goods_data.field04 - HPSdec[i]
   j+= 1
endfor

gdist = dblarr(144, n_elements(goods_data.field03))
for i = 0,143 do begin
   gdist[i, *] = sqrt(cos(gval2[i,*] * !pi/180d)^2 * gval1[i, *]^2 + gval2[i, *]^2)
endfor

goods_min = intarr(144)
for i = 0, 143 do begin
   goods_min[i] = where(gdist[i,*] eq min(gdist[i, *]))   ; returns indices of cos_data array of the likely objects
endfor

openpps, 'z_comp/corresp_pts_GOODSN_zspec'
plot, goods_data.field03[goods_min], goods_data.field04[goods_min], psym=6, yrange=[62.17, 62.26], ystyle=1, title='GOODS-N (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)'
oplot, HPSRA[324:467], HPSdec[324:467], psym=4, color=240
closepps

goods_z = goods_data.field36[goods_min]
new2_z =  z[324:467]
yay = where(goods_z gt 0 and goods_z lt 5, complement = boo)
remove, boo, goods_z, new2_z

openpps, 'z_comp/HPS_GOODSN_zspec'
plot, goods_data.field36[goods_min], z[324:467], psym=6, xrange=[0, 3.5], xstyle=1, yrange=[0, 3.5], ystyle=1, ytitle='HPS redshifts', xtitle='GOODS-N spec z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(goods_z, new2_z)))
closepps

survey_z = [cosmos_z, goods_z]
new_z =  [new1_z, new2_z]

openpps, 'z_comp/HPS_COSGOODS_zspec'
plot, survey_z, new_z, psym=6, xrange=[0, 3.6], yrange=[0, 3.6], ystyle=1, xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS + GOODS-N photo z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(survey_z, new_z)))
closepps
stop

end
