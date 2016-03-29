pro morph_SFR_mass

; Note: These don't include the galaxies that had empty
; outmorph.dat folders, those have been removed
readcol, 'morphology/G_M20.txt', dir, G, M20, R_circ, flag, format = 'I,F,F,F,F'
readcol, 'morphology/mass_SFR.txt', ID, mass, SFR, EBV, z, format = 'I,D,D, D,D'

; Cutout the X-ray bright sources
AGN =  [239, 267, 366, 369, 410, 429, 464] ; AGNs from COSMOS+GOODS
ind = intarr(n_elements(AGN))
for i = 0, n_elements(AGN)-1 do begin
   ind[i] = where(dir eq AGN[i])
endfor

remove, ind, dir, G, M20, flag, R_circ

; Remove stupid galaxy that didn't make it through the R-1 cutoff
ind = where(dir eq 386)
remove, ind, dir, G, M20, flag, R_circ

; Match up what made it through the morphology with what made it
; through the SED fitting
ind = intarr(n_elements(dir))
for i = 0, n_elements(dir)-1 do begin
   ind[i] = where(ID eq dir[i])
endfor

ID = ID[ind]
mass = mass[ind]
SFR = SFR[ind]
z = z[ind]

; Now have to remove flag = 1 galaxies
ind = where(flag eq 1)

remove, ind, G, M20, mass, SFR, R_circ, z, ID

print, 'Pearson correlation coefficient (mass, G):', correlate(mass, G)
print, 'Pearson correlation coefficient (mass, M20):', correlate(mass, M20)
print, 'Pearson correlation coefficient (SFR, G):', correlate(SFR, G)
print, 'Pearson correlation coefficient (SFR, M20):', correlate(SFR, M20)
print, 'Pearson correlation coefficient (SFR, M20):', correlate(SFR-mass, M20)
print, 'Pearson correlation coefficient (SFR, G):', correlate(SFR-mass, G)

b = dblarr(n_elements(G))
for  i = 0, n_elements(G)-1 do begin
   if (G[i] gt 0.45) and (M20[i] gt -2) then begin
       b[i] = dir[i]
    endif
endfor
;print, b

loadct, 17

;plt.vlines(-2.5, 0.18, 0.25)
;plt.hlines(0.215, -2.6, -2.4)

openpps, 'SFR_plots//HPS_o2sed_IRAC/MassvG'
;!P.multi = [0, 3, 2]
multiplot, [3,2]
plot, mass, G, ytitle = 'G', psym = 3, xtickformat = '(A1)', ystyle = 1, yrange = [0.3, 0.8]
oplot, mass, G, psym = sym(1), color = 100, symsize = .8
oploterror, [7.6], [0.73], [0.0697], [0], /lobar      
oploterror, [7.6], [0.73], [0.101], [0], /hibar 
oploterror, [7.6], [0.73], [0], [0.035], /hibar
oploterror, [7.6], [0.73], [0], [0.035], /lobar
multiplot
plot, SFR, G, ytickformat = '(A1)', psym = 3, xtickformat = '(A1)', ystyle = 1, yrange = [0.3, 0.8], xstyle = 1, xrange = [-2, 2.2]
oplot, SFR, G, psym = sym(1), color = 100, symsize = .8
oploterror, [-1.4], [0.73], [0.105], [0], /lobar      
oploterror, [-1.4], [0.73], [0.088], [0], /hibar 
oploterror, [-1.4], [0.73], [0], [0.035], /hibar
oploterror, [-1.4], [0.73], [0], [0.035], /lobar
multiplot
plot, SFR-mass, G, ytickformat = '(A1)', psym = 3, xtickformat = '(A1)', ystyle = 1, yrange = [0.3, 0.8], xstyle = 1, xrange = [-11, -6.5]
oplot, SFR-mass, G, psym = sym(1), color = 100, symsize = .8
oploterror, [-10.4], [0.73], [0.324], [0], /lobar      
oploterror, [-10.4], [0.73], [0.319], [0], /hibar 
oploterror, [-10.4], [0.73], [0], [0.035], /hibar
oploterror, [-10.4], [0.73], [0], [0.035], /lobar
multiplot
plot, mass, M20, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('M_{20}'), psym = 3, ystyle = 1, yrange = [ -2.75, -.25], xtickname = ['7', '8', '9', '10', '11', ' ']
oplot, mass, M20, psym = sym(1), color = 55, symsize = .8
oploterror, [7.6], [-0.6], [0.0697], [0], /lobar      
oploterror, [7.6], [-0.6], [0.101], [0], /hibar 
oploterror, [7.6], [-0.6], [0], [0.1], /hibar
oploterror, [7.6], [-0.6], [0], [0.1], /lobar
multiplot
plot, SFR, M20, xtitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ytickformat = '(A1)', psym = 3, xstyle = 1, xrange = [-2, 3], ystyle = 1, yrange = [-2.75, -.25], xtickname = [' ', '-1', '0', '1', '2', ' ']
oplot, SFR, M20, psym = sym(1), color = 55, symsize = .8
oploterror, [-1.4], [-0.6], [0.105], [0], /lobar      
oploterror, [-1.4], [-0.6], [0.088], [0], /hibar 
oploterror, [-1.4], [-0.6], [0], [0.1], /hibar
oploterror, [-1.4], [-0.6], [0], [0.1], /lobar
multiplot
plot, SFR-mass, M20, xtitle = textoidl('log(sSFR [yr^{-1}])'), ytickformat = '(A1)', psym = 3, xstyle = 1, xrange = [-11, -6.5], ystyle = 1, yrange = [ -2.75, -.25], xtickname = [' ', '-10', '-9-', '-8', '-7']
oplot, SFR-mass, M20, psym = sym(1), color = 55, symsize = .8
oploterror, [-10.4], [-0.6], [0.324], [0], /lobar      
oploterror, [-10.4], [-0.6], [0.319], [0], /hibar 
oploterror, [-10.4], [-0.6], [0], [0.1], /hibar
oploterror, [-10.4], [-0.6], [0], [0.1], /lobar
multiplot
closepps


sigma = 0.003
mean = 1
array = RANDOMN(seed, (size(G))[1])
noise1 = array * sigma + mean
G = G*noise1
sigma = 0.003
mean = 1
array = RANDOMN(seed, (size(M20))[1])
noise2 = array * sigma + mean
M20 = M20 * noise2


openpps, 'SFR_plots/HPS_o2sed_IRAC/GM20'
plot, [-0.5, -2.75], [0.43, 0.655], xtitle = textoidl('M_{20}'), ytitle = 'G', xstyle = 1, xrange = [-0.25, -3], ystyle = 1, yrange = [0.37, 0.7], linestyle = 3, xtickname=[ '-2.5','-2', '-1.5','-1.0','-0.5','0.0']
oplot, M20, G, psym = sym(1), color = 140
oploterror, [-0.5], [0.65], [0.1], [0.035]
closepps



;plt.scatter(M20, G, color = 'green', linewidth = 1.4)
;ax.xaxis.set_minor_locator(minorLocator)
;plt.plot([-0.5, -2.75], [0.43, 0.655], 'k--', linewidth = 1.4)
;plt.vlines(-0.35, 0.685, 0.615, linewidth = 1.4)
;plt.hlines(0.65, -0.25, -0.45, linewidth = 1.4)
;plt.xlabel('M$_{20}$')
;plt.ylabel('G')
;plt.ylim(0.37, 0.7)
;plt.xlim(-3,0)
;plt.gca().invert_xaxis()
;plt.xticks([0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0], ['0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0'])
;plt.gca().tick_params(width = 1.3)







d = dblarr(n_elements(z))
size = dblarr(n_elements(z))
for i = 0, n_elements(z)-1 do begin
   d[i] = lumdist(z[i], /silent)/(1+z[i])^2
   size[i] = d[i]*R_circ[i]/206.265
endfor

openpps, 'SFR_plots/HPS_o2sed_IRAC/size_hist'
cgHistoplot, size, bin = 0.75, ytitle = 'Number of [O II] Emitters', xtitle = textoidl('Physical size (kpc)'), thick = 6,  /line_fill, color = 'forest green', polycolor =  'forest green', orientation = [45, -45], /oprobability, probcolor = 'black'
;axis, yaxis = 1, ytitle = ' ', ycharsize = '1', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks = 5
closepps

R_circ_old = R_circ
size_old = size
outliers = [139, 142]
;print, R_circ[outliers]
;print, size[outliers]
;print, z[outliers]
;print, ID[outliers]
remove, outliers, size, z, R_circ, mass, ID

;print, percentiles(R_circ, value = [0.16, 0.5, 0.84])
;print, percentiles(size, value = [0.16, 0.5, 0.84])

; Little section to put non evolving lines on plot
red = findgen(600)/1000
old_low = 1.8 ;kpc
old_hi = 5.8
d = lumdist(red, /silent)/(1+red)^2
new_low  = old_low*206.265/d
new_hi = old_hi*206.265/d

; Median z = 0.379
e = 0.1 ; kpc
d_e = lumdist(0.379, /silent)/(1+0.379)^2
new_e = 0.1*206.265/d_e
print, new_e

bit1 = size[where(z lt 0.2)]
bit2 = size[where(z gt 0.5 and z lt 0.6)]
print, 'Median, z < 0.2:', median(bit1, /even)
print, 'Error on median, z < 0.2:', 1.253*stdev(bit1)/sqrt(n_elements(bit1))
print, 'Median, 0.5 < z < 0.6:', median(bit2,/even)
print, 'Error on median, 0.5 < z < 0.6:', 1.253*stdev(bit2)/sqrt(n_elements(bit2))

bit1 = R_circ[where(z lt 0.2)]
bit2 = R_circ[where(z gt 0.2 and z lt 0.3)]
bit3 = R_circ[where(z gt 0.3 and z lt 0.4)]
bit4 = R_circ[where(z gt 0.4 and z lt 0.5)]
bit5 = R_circ[where(z gt 0.5 and z lt 0.6)]
z_ind = [0.15, 0.25, 0.35, 0.45, 0.55]
bit1 = percentiles(bit1, value = [.16, .5, .84])
bit2 = percentiles(bit2, value = [.16, .5, .84])
bit3 = percentiles(bit3, value = [.16, .5, .84])
bit4 = percentiles(bit4, value = [.16, .5, .84])
bit5 = percentiles(bit5, value = [.16, .5, .84])
low = [bit1[0], bit2[0], bit3[0], bit4[0], bit5[0]]
med = [bit1[1], bit2[1], bit3[1], bit4[1], bit5[1]]
hi = [bit1[2], bit2[2], bit3[2], bit4[2], bit5[2]]

openpps, 'SFR_plots/HPS_o2sed_IRAC/sizevz', /portrait
multiplot, [1, 2]
plot, z, size, psym = 3, ytitle = textoidl('Physical size (kpc)')
oplot, z, size, psym = sym(1), symsize =1, color = 140
multiplot
plot, z, R_circ, psym = 3, ytitle =  textoidl('Angular size (arcseconds)'), ystyle = 1, yrange = [0, 2.5], xtitle = 'Redshift', ytickname = ['0.0', '0.5', '1.0', '1.5', '2.0', ' ']
oplot, z, R_circ, psym = sym(1), symsize =1, color = 140
oplot, red, new_low
oplot, red, new_hi
oplot, z_ind, low, lines = 3
oplot, z_ind, med, lines = 3
oplot, z_ind, hi, lines = 3
closepps


bit1 = (size)[where(mass lt 8.25)]
bit2 = (size)[where(mass ge 8.25 and mass lt 8.8)]
bit3 = (size)[where(mass ge 8.8 and mass lt 9.35)]
bit4 = (size)[where(mass ge 9.35 and mass lt 9.8)]
bit5 = (size)[where(mass ge 9.8 and mass lt 10.25)]
bit6 = (size)[where(mass ge 10.25 and mass lt 10.8)]
bit7 = (size)[where(mass ge 10.8 and mass lt 11.8)]

bit1 = percentiles(bit1, value = [0.16, 0.5, 0.84])
bit2 = percentiles(bit2, value = [0.16, 0.5, 0.84])
bit3 = percentiles(bit3, value = [0.16, 0.5, 0.84])
bit4 = percentiles(bit4, value = [0.16, 0.5, 0.84])
bit5 = percentiles(bit5, value = [0.16, 0.5, 0.84])
bit6 = percentiles(bit6, value = [0.16, 0.5, 0.84])
bit7 = percentiles(bit7, value = [0.16, 0.5, 0.84])

low = [bit1[0], bit2[0], bit3[0], bit4[0], bit5[0], bit6[0], bit7[0]]
med = [bit1[1], bit2[1], bit3[1], bit4[1], bit5[1], bit6[1], bit7[1]]
high = [bit1[2], bit2[2], bit3[2], bit4[2], bit5[2], bit6[2], bit7[2]]
mass_ind_EW = [7.7, 8.45, 8.95, 9.45, 9.8, 10.35, 10.85]
a = linfit(mass, size)

openpps, 'SFR_plots/HPS_o2sed_IRAC/massvsize'
plot, mass, size, psym = 3,  symsize = 0.9,  ytitle = textoidl('Physical size (kpc)'), xtitle = textoidl('log(M_* [M_{sun}])'), ystyle = 1, yrange = [0.2, 20], xstyle = 1, xrange = [7.5, 11], /ylog
oplot, mass, (size), psym = sym(1),  symsize = 0.9, color = 140
;oplot, [7.5, 11.9], alog10(a[1]*mass+a[0]), line = 4
;oplot, mass_ind_EW, low, lines = 3
;oplot, mass_ind_EW, high, lines = 3
closepps


; For Greg to have sizes                                                                                           
openw, lun, 'morphology/OII_sizes.txt', /get_lun
for i = 0, n_elements(mass)-1 do begin
   printf, lun, string(ID[i])+'   '+string(G[i],format='(F0.2)')+'   '+string(M20[i],format='(F0.2)')+'  '+string(size[i],format='(F0.3)')
endfor
Free_lun, lun


stop
end
