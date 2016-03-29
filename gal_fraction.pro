;+
;
; This program finds the [OII] galaxies in the HPS survey and compares
; that number to the total galaxies at z > 0.5 from COSMOS and GOODS-N
;
;-

pro gal_fraction

readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec

; This section does the COSMOS survey

restore, 'ascii_templates/cosmos_template.sav'
cos_data1 = read_ascii('data/cosmos_zphot_chunk1.tbl', template = cosmos_template)
cos_data2 = read_ascii('data/cosmos_zphot_chunk2.tbl', template = cosmos_template)
cos_data3 = read_ascii('data/cosmos_zphot_chunk3.tbl', template = cosmos_template)

cosRA = [cos_data1.field03, cos_data2.field03, cos_data3.field03]
cosdec = [cos_data1.field04, cos_data2.field04, cos_data3.field04]
cosz = [cos_data1.field05, cos_data2.field05, cos_data3.field05]

loadct, 39
;window, 0
;plot, HPSRA[140:323], HPSdec[140:323], psym = 6, thick=3, yrange=[2.15, 2.4], ystyle=1, xrange=[150.0, 150.22], xstyle=1, title = 'COSMOS (red) and HPS (white) (COSMOS Archive)', xtitle = 'RA (deg)', ytitle = 'dec (deg)'
;oplot, cosRA, cosdec, psym=4, color=240

index = where(cosz gt 0 and cosz le 0.56)
cos_z = cosz[index]
cos_RA = cosRA[index]
cos_dec = cosdec[index]

;window, 1
;plot, HPSRA[140:323], HPSdec[140:323], psym = 6, thick=3, yrange=[2.15, 2.4], ystyle=1, xrange=[150.0, 150.22], xstyle=1, xtitle = 'RA (deg)', ytitle = 'dec (deg)'
;oplot, cos_RA, cos_dec, psym=4, color=240

restore, 'ascii_templates/HPS4_template.sav'
HPS_data = read_ascii('data/HPStable4.txt', template=template)

oxy = HPS_data.field11
remove, ind, oxy

O_index = where(oxy eq '[OII]')
cos_OII = O_index[where(O_index ge 140 and O_index le 323)]
RA = HPSRA[140:323]
dec = HPSdec[140:323]
OIIRA = RA[cos_OII - 140]
OIIdec = dec[cos_OII - 140]

window, 2
plot, OIIRA, OIIdec, psym = 6, thick=3, yrange=[2.19, 2.36], ystyle=1, xrange=[150.02, 150.205], xstyle=1, title = 'COSMOS z < 0.5 (red) and HPS [OII] emitters (white)', xtitle = 'RA (deg)', ytitle = 'dec (deg)'
oplot, cos_RA, cos_dec, psym=4, color=240

print, '-- COSMOS Survey --'
print, 'Number of [OII] galaxies:', n_elements(OIIRA)
print, 'Number of total galaxies:', n_elements(cos_RA)
print, 'Fraction of galaxies that are star-forming:', (float(n_elements(OIIRA))/float(n_elements(cos_RA)))


; This section does the GOODS-N survey

restore, 'ascii_templates/goods_temp.sav'
goods_data = read_ascii('data/goodsn_multiwave_v2.0.cat', template = goods_template)

b = where(goods_data.field03 gt 189.059 and goods_data.field03 lt 189.3065 and goods_data.field04 gt 62.169 and goods_data.field04 lt 62.25503)

;window, 3
;plot, goods_data.field03[b], goods_data.field04[b], psym=6, xrange = [189.05, 189.32], xstyle=1, ystyle=1, yrange=[62.165, 62.26], title = 'GOODS-N (white) and HPS (red)', xtitle='RA (deg)', ytitle='dec (deg)'
;oplot, HPSRA[324:467], HPSdec[324:467], psym=4, color=240, thick=3 

goodsRA = goods_data.field03[b]
goodsdec = goods_data.field04[b]

goodsz = goods_data.field38[b]
index = where(goodsz ge 0 and goodsz le 0.56)
goods_z =goodsz[index]
goods_RA = goodsRA[index]
goods_dec = goodsdec[index]

goods_OII = O_index[where(O_index ge 324)]
RA = HPSRA[324:467]
dec = HPSdec[324:467]
OIIRA = RA[goods_OII - 323]
OIIdec = dec[goods_OII - 323]

window, 4
plot, OIIRA, OIIdec, psym = 6, thick=3, xrange = [189.05, 189.32], xstyle=1, ystyle=1, yrange=[62.165, 62.26], title = 'GOODS-N z < 0.5 (red) and HPS [OII] emitters (white)', xtitle = 'RA (deg)', ytitle = 'dec (deg)'
oplot, goods_RA, goods_dec, psym=4, color=240

print, '-- GOODS-N Survey --'
print, 'Number of [OII] galaxies:', n_elements(OIIRA)
print, 'Number of total galaxies:', n_elements(goods_RA)
print, 'Fraction of galaxies that are star-forming:', (float(n_elements(OIIRA))/float(n_elements(goods_RA)))

stop
end
