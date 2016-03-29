;+
;
; Compares spec z's of COSMOS and GOODS-N to HPS redshifts
; (except only some specs z's of COSMOS are available, so the
; rest are photo z's)
;
;-
pro comp_zspec_2

restore, 'ascii_templates/HPS4_template.sav'
data = read_ascii('data/HPStable4.txt', template=template)
z = data.field12

; This section is to acquire what spectroscopic redshifts are
; available for COSMOS

cat = mrdfits('data/ACS-GC_published_catalogs/cosmos_i_public_catalog_V1.0.fits.gz',1)
readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec, z

newRA = cat.RA
newdec = cat.dec
newz = cat.specz

val1_spec = dblarr(184, n_elements(cat.RA))
val2_spec = dblarr(184, n_elements(cat.dec))
j=0
for i = 140, 323 do begin
   val1_spec[j, *] = newRA - HPSRA[i]
   val2_spec[j, *] = newdec - HPSdec[i]
   j+= 1
endfor

dist_spec = dblarr(184, n_elements(newRA))
for i = 0,183 do begin
   dist_spec[i, *] = sqrt(cos(val2_spec[i,*] * !pi/180d)^2 * val1_spec[i, *]^2 + val2_spec[i, *]^2)
endfor

cos_min_spec = lonarr(184)
min_spec = dblarr(184)
for i = 0, 183 do begin
   cos_min_spec[i] = where(dist_spec[i,*] eq min(dist_spec[i, *]))   ; returns indices of cos_data array of the likely objects
   min_spec[i] = min(dist_spec[i, *])
endfor

loadct, 39
openpps, 'z_comp/corresp_pts_COSMOS_zspec_ACS'   
plot, newRA[cos_min_spec], newdec[cos_min_spec], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
oplot, HPSRA[140:323], HPSdec[140:323], psym=4, color=240
closepps

ACS_z = newz[cos_min_spec]
new_z =  z[140:323]
yay = where(ACS_z gt 0 and ACS_z lt 5, complement = bad)
remove, bad, ACS_z, new_z 

openpps, 'z_comp/HPS_COSMOS_zspec_ACS'
plot, ACS_z, new_z, psym=6, xrange=[0, 4], xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS spec z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(ACS_z, new_z)))
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
new_z =  z[324:467]
yay = where(goods_z gt 0 and goods_z lt 5, complement = boo)
remove, boo, goods_z, new_z

openpps, 'z_comp/HPS_GOODSN_zspec'
plot, goods_data.field36[goods_min], z[324:467], psym=6, xrange=[0, 3.5], xstyle=1, yrange=[0, 3.5], ystyle=1, ytitle='HPS redshifts', xtitle='GOODS-N spec z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(goods_z, new_z)))
closepps

survey_z = [newz[cos_min_spec], goods_data.field36[goods_min]]
new_z =  z[140:467]
yay = where(survey_z gt 0 and survey_z lt 5, complement = bad)
remove, bad, survey_z, new_z 

openpps, 'z_comp/HPS_COSGOODS_zspec_ACS'
plot, survey_z, new_z, psym=6, xrange=[0, 0.6], yrange=[0, .6], ystyle=1, xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS + GOODS-N spec+photo z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(survey_z, new_z)))
closepps
stop

end
