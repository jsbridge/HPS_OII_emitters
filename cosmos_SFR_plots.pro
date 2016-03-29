pro cosmos_SFR_plots

cosmos_OII = cosmos_OII()

; Note the use of .likestats files - they have the best fit, as
; opposed to .margestats, which
; have mean fits
restore, 'ascii_templates/likestats_template_tau.sav'   ; will need to change template if not doing exponential SFH
file = findfile('SED_fitting/COSMOS_SEDfits/o2sed_40ktau/*.likestats')  ; modify the folder name for the different fits

mass = dblarr(n_elements(file))
SFR = dblarr(n_elements(file))
age = dblarr(n_elements(file))

; This program calls the SED fits.  I've just been modifying it
; to work for various versions and paramters

for i = 0, n_elements(file)-1 do begin
   data = read_ascii(file[i], template = likestats_template_tau)   ;  Check you've got the right template
   mass[i] = data.field2[4]    ; Check the subscripting numbers here - depending on what params have been fit, they may be in different rows or whatnot
;   SFR[i] = data.field2[3]/data.field2[0]           ; This SFR is for constant SFH
   SFR[i] = data.field2[4]/(2.71828182846^(1d/data.field2[3]))  ; This SFR is for exponential SFH
   age[i] = data.field2[0]
endfor

mass = mass*alog10(2.71828182846)   ; log(mass [M_sol])
SFR = SFR*alog10(2.71828182846)       ; log(SFR [M_sol/yr])
SSFR = alog10(10^(SFR)/10^(mass))     ; log(SSFR [yr^-1])
age = age*alog10(2.71828182846)      ; log(age)



; This following section is for if not all 110 galaxies made it
; through.  It basically pulls on the surviving z's from the
; cosmos_OII structure.  If all galaxies made it through, comment out
; this section

;ind = fix(strmid(file, 46, 3))  ;Change '47' depending on lenght of file names!!!!!
;x = intarr(n_elements(ind))
;for i = 0, n_elements(x)-1 do begin
;   x[i] = where(cosmos_OII.id eq ind[i])
;endfor
;y = indgen(n_elements(cosmos_OII.id))
;remove, x, y
z = cosmos_OII.z
;remove, y, z


; This section splits up the SFR and mass by redshift, for extra nice plotting
SFR_lowz = dblarr(n_elements(SFR))
SFR_mid1z = dblarr(n_elements(SFR))
SFR_mid2z = dblarr(n_elements(SFR))
SFR_hiz = dblarr(n_elements(SFR))
mass_lowz = dblarr(n_elements(mass))
mass_mid1z = dblarr(n_elements(mass))
mass_mid2z = dblarr(n_elements(mass))
mass_hiz =  dblarr(n_elements(mass))
lowz =  dblarr(n_elements(mass))
mid1z =  dblarr(n_elements(mass))
mid2z =  dblarr(n_elements(mass))
hiz =  dblarr(n_elements(mass))
SSFR_lowz = dblarr(n_elements(SSFR))
SSFR_mid1z = dblarr(n_elements(SSFR))
SSFR_mid2z = dblarr(n_elements(SSFR))
SSFR_hiz = dblarr(n_elements(SSFR))

for i = 0, n_elements(z)-1 do begin
   if z[i] ge 0 and z[i] lt 0.2 then begin
      SFR_lowz[i] = SFR[i]
      mass_lowz[i] = mass[i]
      lowz[i] = z[i]
      SSFR_lowz[i] = SSFR[i]
   endif
   if z[i] ge 0.2 and z[i] lt 0.32 then begin
      SFR_mid1z[i] = SFR[i]
      mass_mid1z[i] = mass[i]
      mid1z[i] = z[i]
      SSFR_mid1z[i] = SSFR[i]
   endif
   if z[i] ge 0.32 and z[i] lt 0.45 then begin
      SFR_mid2z[i] = SFR[i]
      mass_mid2z[i] = mass[i]
      mid2z[i] = z[i]
      SSFR_mid2z[i] = SSFR[i]
   endif
   if z[i] ge 0.45 and z[i] le 0.57 then begin
      SFR_hiz[i] = SFR[i]
      mass_hiz[i] = mass[i]
      hiz[i] = z[i]
      SSFR_hiz[i] = SSFR[i]
   endif
endfor

; The following section reads in Robin's line [OII] SFRs (as
; well as UV SFRs)
readcol, 'data/robin_sfr.dat', robin_id, robin_RA, robin_dec, robin_UV, robin_OII, format = 'I, D, D, D, D'

; Expand this to include GOODS!
arr = intarr(n_elements(cosmos_OII.id))
for i = 0, n_elements(cosmos_OII.id) - 1 do begin
   arr[i] = findel(robin_id, cosmos_OII.id[i])
endfor

; Divide mass into chunks, find the median and +/- stddev for each bin
; for the purpose of plotting
bit1 = (robin_OII[arr])[where(mass lt 7.67)]
bit2 = (robin_OII[arr])[where(mass ge 7.67 and mass lt 8.33)]
bit3 = (robin_OII[arr])[where(mass ge 8.33 and mass lt 9)]
bit4 = (robin_OII[arr])[where(mass ge 9 and mass lt 9.67)]
bit5 = (robin_OII[arr])[where(mass ge 9.67 and mass lt 10.33)]
bit6 = (robin_OII[arr])[where(mass ge 10.33 and mass lt 11)]

std1 = stddev(bit1)
std2 = stddev(bit2)
std3 = stddev(bit3)
std4 = stddev(bit4)
std5 = stddev(bit5)
std6 = stddev(bit6)

med1 = median(bit1)
med2 = median(bit2)
med3 = median(bit3)
med4 = median(bit4)
med5 = median(bit5)
med6 = median(bit6)

med = [med1, med2, med3, med4, med5, med6]
std_up = [med1+std1, med2+std2, med3+std3, med4+std4, med5+std5, med6+std6]
std_low = [med1-std1, med2-std2, med3-std3, med4-std4, med5-std5, med6-std6]
mass_ind = [7.33, 8, 8.67, 9.33, 10, 10.67]

other1 = (robin_OII[arr]-mass)[where(mass lt 7.67)]
other2 = (robin_OII[arr]-mass)[where(mass ge 7.67 and mass lt 8.33)]
other3 = (robin_OII[arr]-mass)[where(mass ge 8.33 and mass lt 9)]
other4 = (robin_OII[arr]-mass)[where(mass ge 9 and mass lt 9.67)]
other5 = (robin_OII[arr]-mass)[where(mass ge 9.67 and mass lt 10.33)]
other6 = (robin_OII[arr]-mass)[where(mass ge 10.33 and mass lt 11)]

ostd1 = stddev(other1)
ostd2 = stddev(other2)
ostd3 = stddev(other3)
ostd4 = stddev(other4)
ostd5 = stddev(other5)
ostd6 = stddev(other6)

omed1 = median(other1)
omed2 = median(other2)
omed3 = median(other3)
omed4 = median(other4)
omed5 = median(other5)
omed6 = median(other6)

omed = [omed1, omed2, omed3, omed4, omed5, omed6]
ostd_up = [omed1+ostd1, omed2+ostd2, omed3+ostd3, omed4+ostd4, omed5+ostd5, omed6+ostd6]
ostd_low = [omed1-ostd1, omed2-ostd2, omed3-ostd3, omed4-ostd4, omed5-ostd5, omed6-ostd6]

; Import Pirzkal and Noeska data for comparison
readcol, 'data/PN_Mass_SFR/Kai0.2-0.45.txt', N_mass1, N_SFR1, N_mass_up1, N_SFR_up1, N_mass_low1, N_SFR_low1, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.45-0.70.txt', N_mass2, N_SFR2, N_mass_up2, N_SFR_up2, N_mass_low2, N_SFR_low2, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.70-0.85.txt', N_mass3, N_SFR3, N_mass_up3, N_SFR_up3, N_mass_low3, N_SFR_low3, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.85-1.10.txt', N_mass4, N_SFR4, N_mass_up4, N_SFR_up4, N_mass_low4, N_SFR_low4, format = 'D, D, D, D, D, D'

N_mass_init = [N_mass1, N_mass2];, N_mass3, N_mass4]
N_mass = N_mass_init[sort(N_mass_init)]
N_SFR = [N_SFR1, N_SFR2, N_SFR3, N_SFR4]
N_SFR = N_SFR[sort(N_mass_init)]
; N_mass_up = [N_mass_up1, N_mass_up2, N_mass_up3, N_mass_up3]
N_SFR_up = [N_SFR_up1, N_SFR_up2, N_SFR_up3, N_SFR_up4]
N_SFR_up = N_SFR_up[sort(N_mass_init)]
; N_mass_low = [N_mass_low1, N_mass_low2, N_mass_low3, N_mass_low4]
N_SFR_low = [N_SFR_low1, N_SFR_low2, N_SFR_low3, N_SFR_low4]
N_SFR_low = N_SFR_low[sort(N_mass_init)]

readcol, 'data/PN_Mass_SFR/SFR_Mass_0.2_0.5.txt', P_mass1, P_SFR1, P_err1, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_0.5_0.7.txt', P_mass2, P_SFR2, P_err2, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_0.7_0.8.txt', P_mass3, P_SFR3, P_err3, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_0.8_1.1.txt', P_mass4, P_SFR4, P_err4, format = 'D, D, D'

P_mass_init = [P_mass1, P_mass2, P_mass3, P_mass4]
P_mass = P_mass_init[sort(P_mass_init)]
P_SFR = [P_SFR1, P_SFR2, P_SFR3, P_SFR4]
P_SFR = P_SFR[sort(P_mass_init)]
P_SFR_up = [P_SFR1+P_err1, P_SFR2+P_err2, P_SFR3+P_err3, P_SFR4+P_err4]
P_SFR_up = P_SFR_up[sort(P_mass_init)]
P_SFR_low = [P_SFR1-P_err1, P_SFR2-P_err2, P_SFR3-P_err3, P_SFR4-P_err4]
P_SFR_low = P_SFR_low[sort(P_mass_init)]

; PLOT ALL THE THINGS!!! Note that I have to remember to have moved the
; old plots into a different folder, or they will be written over with
; these
A = FINDGEN(36) *  (!PI*3/36.) 
USERSYM, COS(A), SIN(A), /FILL 
loadct, 39
openpps, 'SFR_plots/cosmos_lineSFRvsmass_z'
plot, mass_lowz, robin_OII[arr], psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7, 11], ystyle = 1, yrange = [-2, 3]
oplot, mass_lowz, robin_OII[arr], psym = 2, color = 185
oplot, mass_mid1z, robin_OII[arr], psym = 6, color = 220
oplot, mass_mid2z, robin_OII[arr], psym = 5, color = 150
oplot, mass_hiz, robin_OII[arr], psym = 4, color = 50
oplot, mass_ind, med, psym = 8, symsize = 1.9
oplot, mass_ind, std_up
oplot, mass_ind, std_low
;oplot, N_mass, N_SFR, psym = 8, symsize = 1.1, color = 30
;oplot, N_mass, N_SFR_up, color = 30, linestyle = 3
;oplot, N_mass, N_SFR_low, color = 30, linestyle = 3
;oplot, P_mass, P_SFR, psym = 8, color =105
;oplot, P_mass, P_SFR_up, color = 105, linestyle = 3
;oplot, P_mass, P_SFR_low, color = 105, linestyle = 3
legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [2, 6, 5, 4], color = [185, 220, 150, 50]
closepps


openpps, 'SFR_plots/cosmos_lineSFRvsmass_z_PN'
!p.multi = [0, 3, 2]
plot, mass_lowz, robin_OII[arr], psym = 5, ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [5, 12], ystyle = 1, yrange = [-2, 3], pos = [.05, .5, .35, 0.95], XTickformat='(A1)', ytickname = [' ', '-1', '0', '1', '2', '3']
xyouts, 9, 2.5,  '0.0 < z < 0.2', charsize = .9
legend, ['This survey'], psym = [5], charsize = .75

plot, mass_mid1z, robin_OII[arr], psym = 5, xstyle = 1, xrange = [5, 12], ystyle = 1, yrange = [-2, 3], pos = [.35, .5, .65, .95], XTickformat='(A1)', YTickformat='(A1)'
oplot, mass_mid2z, robin_OII[arr], psym = 5
oplot, P_mass1, P_SFR1, psym = 6, color = 50
oplot, N_mass1, N_SFR1, psym = 4, color = 220
oplot, N_mass1, N_SFR_low1, lines = 0, color = 220
oplot, N_mass1, N_SFR_up1, lines=0, color = 220
oplot, P_mass1, P_SFR1+P_err1, lines=0, color = 50
oplot, P_mass1, P_SFR1-P_err1, lines=0, color=50
xyouts, 9, 2.5,  '0.2 < z < 0.45', charsize = .9
legend, [ 'Noeske et al. 2007', 'Pirzkal et a. 2012', 'This survey'], psym = [ 4, 6, 5], color = [220, 50, 0], charsize = .75

plot, mass_hiz, robin_OII[arr], psym = 5, xtitle = textoidl('log(M_* [M_{sun}])'), xstyle = 1, xrange = [5, 12], ystyle = 1, yrange = [-2, 3], pos = [.65, .5, .95, .95], YTickformat='(A1)', xtickname=[' ', '6', '7', '8', '9', '10', '11', '12']
oplot, N_mass2, N_SFR2, psym = 4, color = 220
oplot, P_mass2, P_SFR2, psym = 6, color = 50
oplot, N_mass2, N_SFR_low2, lines = 0, color = 220
oplot, N_mass2, N_SFR_up2, lines=0, color = 220
oplot, P_mass2, P_SFR2+P_err2, lines=0, color = 50
oplot, P_mass2, P_SFR2-P_err2, lines=0, color=50
xyouts, 9, 2.5,  '0.45 < z < 0.7', charsize = .9
legend, [ 'Noeske et al. 2007', 'Pirzkal et a. 2012', 'This survey'], psym = [ 4, 6, 5], color = [220, 50, 0], charsize = .75

plot, N_mass3, N_SFR3, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [5, 12], ystyle = 1, yrange = [-2, 3], pos = [.05, .1, .35, .5], xtickname=['5', '6', '7', '8', '9', '10', '11', ' '], ytickname = ['-2', '-1', '0', '1', '2', ' ']
oplot, N_mass3, N_SFR3, psym = 4, color = 220
oplot, P_mass3, P_SFR3, psym = 6, color = 50
xyouts, 9, 2.5,  '0.7 < z < 0.85', charsize = .9
oplot, N_mass3, N_SFR_low3, lines = 0, color = 220
oplot, N_mass3, N_SFR_up3, lines=0, color = 220
oplot, P_mass3, P_SFR3+P_err3, lines=0, color = 50
oplot, P_mass3, P_SFR3-P_err3, lines=0, color=50
legend, ['Noeske et al. 2007', 'Pirzkal et a. 2012'], psym = [ 4, 6], color = [ 220, 50], charsize = .75

plot, N_mass4, N_SFR4, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), xstyle = 1, xrange = [5, 12], ystyle = 1, yrange = [-2, 3], pos = [.35, .1, .65, .5] , YTickformat='(A1)', xtickname=[' ', '6', '7', '8', '9', '10', '11', '12']
oplot, N_mass4, N_SFR4, color = 220, psym = 4
oplot, P_mass4, P_SFR4, psym = 6, color = 50
xyouts, 9, 2.5,  '0.85 < z < 1.1', charsize = .9
oplot, N_mass4, N_SFR_low4, lines = 0, color = 220
oplot, N_mass4, N_SFR_up4, lines=0, color = 220
oplot, P_mass4, P_SFR4+P_err4, lines=0, color = 50
oplot, P_mass4, P_SFR4-P_err4, lines=0, color=50
legend, ['Noeske et al. 2007', 'Pirzkal et a. 2012'], psym = [ 4, 6], color = [ 220, 50], charsize = .75

closepps


openpps, 'SFR_plots/cosmos_lineSSFRvsmass_z'
plot, mass_lowz, robin_OII[arr]-mass_lowz, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SSFR [yr^{-1}])'), xstyle = 1, xrange = [7, 12], ystyle = 1, yrange =[-11, -5]
oplot, mass_lowz, robin_OII[arr]-mass_lowz, psym = 2, color = 185
oplot, mass_mid1z, robin_OII[arr]-mass_mid1z, psym = 6, color = 220
oplot, mass_mid2z, robin_OII[arr]-mass_mid2z, psym = 5, color = 150
oplot, mass_hiz, robin_OII[arr]-mass_hiz, psym = 4, color = 50
oplot, mass_ind, omed, psym = 8, symsize = 1.9
oplot, mass_ind, ostd_up
oplot, mass_ind, ostd_low
legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [2, 6, 5, 4], color = [185, 220, 150, 50]
closepps

openpps, 'SFR_plots/cosmos_mass_hist'
cgHistoplot, mass, xtitle = textoidl('log(M_* [M_{sol}])'), ytitle = 'Number of [OII] Galaxies in COSMOS', bin = 0.2, /line_fill, color = 'red', thick = 3, polycolor = 'red', orientation = [45, -45]
closepps

openpps, 'SFR_plots/cosmos_lineSRF_hist'
cgHistoplot, robin_OII[arr], xtitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ytitle = 'Number of [OII] Galaxies in COSMOS', bin = 0.2, /line_fill, color = 'blue', thick = 3, polycolor = 'blue', orientation = [45, -45]
closepps

openpps, 'SFR_plots/cosmos_lineSSRF_hist'
cgHistoplot, robin_OII[arr]-mass, xtitle = textoidl('log(SSFR [yr^{-1}])'), ytitle = 'Number of [OII] Galaxies in COSMOS', bin = 0.2, /line_fill, color = 'orange', thick = 3, polycolor = 'orange', orientation = [45, -45]
closepps

readcol, 'data/robin_completeness.dat', complete_z, complete_90, complete_80, complete_70, complete_60, complete_50, format = 'D,D,D,D,D,D'
complete_90_SFR = complete_90 - 41.1817d
complete_80_SFR = complete_80 - 41.1817d
complete_70_SFR = complete_70 - 41.1817d
complete_60_SFR = complete_60 - 41.1817d
complete_50_SFR = complete_50 - 41.1817d

openpps, 'SFR_plots/cosmos_lineSFRvsz' 
plot,z, (robin_OII[arr]), psym = 3, xtitle = textoidl('Redshift (z)'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-2, 2.5], xstyle = 1, xrange = [0.02, 0.58]
oplot, z, robin_OII[arr], psym = 6, color = 205
oplot, complete_z, complete_80_SFR
closepps




; The following plots are all geared up to use the SED SFR, which
; currently doesn't work, so it's all commented out

;openpps, 'SFR_plots/cosmos_SFRvsmass_z'
;plot, mass_lowz, SFR_lowz, psym = 2, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7, 12], ystyle = 1, yrange = [-.1, 0.8]
;oplot, mass_mid1z, SFR_mid1z, psym = 6, color = 220
;oplot, mass_mid2z, SFR_mid2z, psym = 5, color = 150
;oplot, mass_hiz, SFR_hiz, psym = 4, color = 50
;legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [2, 6, 5, 4], color = [0, 220, 150, 50]
;closepps

;openpps, 'SFR_plots/cosmos_SSFRvsmass_z'
;plot, mass_lowz, SSFR_lowz, psym = 2, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SSFR [yr^{-1}])'), xstyle = 1, xrange = [7, 12], ystyle = 1, yrange =[-11, -5]
;oplot, mass_mid1z, SSFR_mid1z, psym = 6, color = 220
;oplot, mass_mid2z, SSFR_mid2z, psym = 5, color = 150
;oplot, mass_hiz, SSFR_hiz, psym = 4, color = 50
;egend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [2, 6, 5, 4], color = [0, 220, 150, 50]
;closepps

;openpps, 'SFR_plots/cosmos_SFRvsz' 
;plot, z, SFR, psym = 6, xtitle = textoidl('Redshift (z)'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])');, ystyle = 1, yrange = [0.2, 0.8]
;closepps



stop
end
