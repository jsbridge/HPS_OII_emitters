pro HPSfield

data = read_ascii('HPStable3.txt', comment_symbol=';')

RA = 15d * (data.field01[1,*] + data.field01[2,*]/60d + data.field01[3,*]/3600d)
pos_dec = data.field01[4,41:478] + data.field01[5,41:478]/60d + data.field01[6,41:478]/3600d
neg_dec = data.field01[4, 0:40] - data.field01[5,0:40]/60d - data.field01[6,0:40]/3600d
dec = [[neg_dec], [pos_dec]]

;window, 0
;plot, RA[0:40], dec[0:40], psym=6
;window, 1
;plot, RA[41:141], dec[41:141], psym=6
;window, 2
;plot, RA[142:330], dec[142:330], psym=6, yrange=[2.2, 2.35], ystyle=1
;window, 3
;plot, RA[331:478], dec[331:478], psym=6, yrange=[62.16, 62.26],
;ystyle=1

restore, 'allcos_template.sav'
allcos_data = read_ascii('cosmos_zphot_mag25.tbl', template=allcos_template)
loadct, 39
plot, allcos_data.field03, allcos_data.field04, xrange = [149.41, 150.83], xstyle=1, yrange=[1.498, 2.913], ystyle=1, title="COSMOS (white), HPS (red)", ytitle='Dec (deg)', xtitle='RA (deg)'
oplot, RA, dec, psym=6, thick=2, color=255 

print, 'Max RA(1):', max(RA[0:40])
print, 'Min RA(1):', min(RA[0:40])
print, 'Max Dec(1):', max(dec[0:40])
print, 'Min Dec(1):', min(dec[0:40])
print, 'Max RA(2):', max(RA[41:141])
print, 'Min RA(2):', min(RA[41:141])
print, 'Max Dec(2):', max(dec[41:141])
print, 'Min Dec(2):', min(dec[41:141])
print, 'Max RA(3):', max(RA[142:330])
print, 'Min RA(3):', min(RA[142:330])
print, 'Max Dec(3):', max(dec[142:330])
print, 'Min Dec(3):', min(dec[142:330])
print, 'Max RA(4):', max(RA[331:478])
print, 'Min RA(4):', min(RA[331:478])
print, 'Max Dec(4):', max(dec[331:478])
print, 'Min Dec(4):', min(dec[331:478])

stop

end
