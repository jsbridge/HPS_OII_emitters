;+
;
; Compares the z's of HPS to photo z's of (hopefully)
; corresponding opjects in COSMOS and GOODS-N
;
;-

pro comp_zphot

restore, 'ascii_templates/HPS4_template.sav'
data = read_ascii('data/HPStable4.txt', template=template)
z = data.field12

readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec, z

restore, 'ascii_templates/cosmos_template.sav'
cos_data = read_ascii('data/cosmos_zphot_HPSfield.tbl', template = cosmos_template)

val1 = dblarr(184, n_elements(cos_data.field03))
val2 = dblarr(184, n_elements(cos_data.field04))
j=0
for i = 140, 323 do begin
   val1[j, *] = cos_data.field03 - HPSRA[i]
   val2[j, *]= cos_data.field04 - HPSdec[i]
   j+= 1
endfor

dist = dblarr(184, n_elements(cos_data.field03))
for i = 0,183 do begin
   dist[i, *] = sqrt(cos(val2[i,*] * !pi/180d)^2 * val1[i, *]^2 + val2[i, *]^2)
endfor

cos_min = intarr(184)
for i = 0, 183 do begin
   cos_min[i] = where(dist[i,*] eq min(dist[i, *]))   ; returns indices of cos_data array of the likely objects
endfor

loadct, 39
openpps, 'z_comp/corresp_pts_COSMOS_zphot'   
plot, cos_data.field03[cos_min], cos_data.field04[cos_min], psym=6, title='COSMOS (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)', yrange=[2.20, 2.35], ystyle=1, xrange=[150.02, 150.2], xstyle=1
oplot, HPSRA[140:323], HPSdec[140:323], psym=4, color=240
closepps

survey_z = [cos_data.field05[cos_min]]
new_z =  z[140:323]
yay = where(survey_z gt 0 and survey_z lt 5, complement = bad)
remove, bad, survey_z, new_z 

openpps, 'z_comp/HPS_COSMOS_zphot'
plot, cos_data.field05[cos_min], z[140:323], psym=6, xrange=[0, 4], xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS photo z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(survey_z, new_z)))
closepps

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

openpps, 'z_comp/corresp_pts_GOODSN_zphot'
plot, goods_data.field03[goods_min], goods_data.field04[goods_min], psym=6, yrange=[62.17, 62.26], ystyle=1, title='GOODS-N (white) and corresponding HPS (red)', ytitle='Declination (deg)', xtitle='RA (deg)'
oplot, HPSRA[324:467], HPSdec[324:467], psym=4, color=240
closepps

openpps, 'z_comp/HPS_GOODSN_zphot'
plot, goods_data.field38[goods_min], z[324:467], psym=6, xrange=[0, 3.5], xstyle=1, yrange=[0, 3.5], ystyle=1, ytitle='HPS redshifts', xtitle='GOODS-N photo z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(goods_data.field38[goods_min], z[324:467])))
closepps

survey_z = [cos_data.field05[cos_min], goods_data.field38[goods_min]]
new_z =  z[140:467]
yay = where(survey_z gt 0 and survey_z lt 5, complement = bad)
remove, bad, survey_z, new_z 

openpps, 'z_comp/HPS_COSGOODS_zphot'
plot, survey_z, new_z, psym=6, xrange=[0, 3.6], yrange=[0, 3.6], ystyle=1, xstyle=1, ytitle='HPS redshifts', xtitle='COSMOS + GOODS-N photo z'
xyouts, 1.5, 0.5, strcompress('r = '+string(correlate(survey_z, new_z)))
closepps

stop

end
