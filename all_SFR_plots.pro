pro all_SFR_plots

all_OII = all_OII()
R_cutoff = R_cutoff()

restore, 'ascii_templates/likestats_template.sav'   ; will need to change template depending on if doing exponential SFH or constant
;cos_file = findfile('SED_fitting/cos_80k_like/*.likestats')  ; modify the folder name for the different fits
cos_file = findfile('SED_fitting/cosmos_IRAC_like/*.likestats')
;goods_file = findfile('SED_fitting/goods_80k_like/*.likestats')
goods_file = findfile('SED_fitting/goods_final_like/*.likestats')
munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
;xmm_file = findfile('SED_fitting/xmm_80k_like/*.likestats')
xmm_file = findfile('SED_fitting/xmm_IRAC_like/*.likestats')

xmm_file = [xmm_file[0], xmm_file[7], xmm_file[-3:-1], xmm_file[1:6], xmm_file[8:-4]]
munics_file = [munics_file[-35:-1], munics_file[0:24]]

c_mass = dblarr(n_elements(cos_file))
g_mass = dblarr(n_elements(goods_file))
m_mass = dblarr(n_elements(munics_file))
x_mass = dblarr(n_elements(xmm_file))
c_mass_up = dblarr(n_elements(cos_file))
g_mass_up = dblarr(n_elements(goods_file))
m_mass_up = dblarr(n_elements(munics_file))
x_mass_up = dblarr(n_elements(xmm_file))
c_mass_low = dblarr(n_elements(cos_file))
g_mass_low = dblarr(n_elements(goods_file))
m_mass_low = dblarr(n_elements(munics_file))
x_mass_low = dblarr(n_elements(xmm_file))
c_age = dblarr(n_elements(cos_file))
g_age = dblarr(n_elements(goods_file))
m_age = dblarr(n_elements(munics_file))
x_age = dblarr(n_elements(xmm_file))

for i = 0, n_elements(cos_file)-1 do begin
   c_data = read_ascii(cos_file[i], template = likestats_template) 
   c_mass[i] = c_data.field2[3] 
   c_mass_up[i] = c_data.field4[3]-c_data.field2[3]  ; 1 sigma
   c_mass_low[i] = c_data.field2[3]-c_data.field3[3]  ; 1 sigma
   c_age[i] = c_data.field2[0]
endfor
for i = 0, n_elements(goods_file)-1 do begin
   g_data = read_ascii(goods_file[i], template = likestats_template) 
   g_mass[i] = g_data.field2[3]  
   g_mass_up[i] = g_data.field4[3]-g_data.field2[3]  ; 1 sigma
   g_mass_low[i] = g_data.field2[3]-g_data.field3[3]  ; 1 sigma
   g_age[i] = g_data.field2[0]
endfor
for i = 0, n_elements(munics_file)-1 do begin
   m_data = read_ascii(munics_file[i], template = likestats_template) 
   m_mass[i] = m_data.field2[3]  
   m_mass_up[i] = m_data.field4[3]-m_data.field2[3]  ; 1 sigma
   m_mass_low[i] = m_data.field2[3]-m_data.field3[3]  ; 1 sigma
   m_age[i] = m_data.field2[0]
endfor
for i = 0, n_elements(xmm_file)-1 do begin
   x_data = read_ascii(xmm_file[i], template = likestats_template) 
   x_mass[i] = x_data.field2[3]  
   x_mass_up[i] = x_data.field4[3]-x_data.field2[3]  ; 1 sigma
   x_mass_low[i] = x_data.field2[3]-x_data.field3[3]  ; 1 sigma
   x_age[i] = x_data.field2[0]
endfor

mass_old = [x_mass, m_mass, c_mass, g_mass]
mass_old = mass_old*alog10(2.71828182846)   ; log(mass [M_sol])
z_old = [all_OII.xmm.z, all_OII.munics.z, all_OII.cosmos.z,all_OII.goods.z]
IDs_old = [all_OII.xmm.id, all_OII.munics.id, all_OII.cosmos.id, all_OII.goods.id]
EW_old = [all_OII.xmm.EW, all_OII.munics.EW, all_OII.cosmos.EW, all_OII.goods.EW]
EW_old_low = [all_OII.xmm.EW_low, all_OII.munics.EW_low, all_OII.cosmos.EW_low, all_OII.goods.EW_low]
EW_old_up = [all_OII.xmm.EW_up, all_OII.munics.EW_up, all_OII.cosmos.EW_up, all_OII.goods.EW_up]
age_old = [x_age, m_age, c_age, g_age]
age_old = age_old*alog10(2.71828182846)
Rmag =  [all_OII.xmm.Rmag, all_OII.munics.Rmag, all_OII.cosmos.Rmag, all_OII.goods.Rmag]
; error bar arrays from mass SED errors
xerr_up = [x_mass_up, m_mass_up, c_mass_up, g_mass_up]*alog10(2.71828182846)
xerr_low = [x_mass_low, m_mass_low, c_mass_low, g_mass_low]*alog10(2.71828182846)

; The following section reads in Robin's line [OII] SFRs (as
; well as UV SFRs)
readcol, 'data/robin_sfr.dat', ID_robin, E_BV_robin, UVSFR_robin, lineSFR_robin, format = 'I,X,X,D,D,X,D,X,X'

; Read in my SFRs with my SED E(B-V)
lineSFR_old = luminosity()
;lineSFR_old = flux_SFR()

; Read in the chi squared values for the SED fits
chi = chi_squared()

; This section cuts out the X-ray emitters and stores them elsewhere
; It's commented out because I only want to remove the AGN
; candidates, not ALL X-ray sources - see below
;x = all_OII.xray_ind
;xray_mass = mass_old[x]
;xray_lineSFR = lineSFR_old.SFR[x]
;xray_IDs = IDs_old[x]
;line_SFR_old = lineSFR_old.SFR
;remove, x, mass_old, line_SFR_old, IDs_old, z_old

; R-1 cutoff from SED fits, currently set at 0.2
lum = dblarr(n_elements(lineSFR_old.lum))
mass = dblarr(n_elements(mass_old))
IDs = intarr(n_elements(IDs_old))
z = dblarr(n_elements(z_old))
lineSFR = dblarr(n_elements(lineSFR_old.SFR))
lineSFR_up = dblarr(n_elements(lineSFR_old.SFR_up))
lineSFR_low = dblarr(n_elements(lineSFR_old.SFR_low))
SFR_red = dblarr(n_elements(lineSFR_old.SFR_red))
EW = dblarr(n_elements(EW_old))
EW_up = dblarr(n_elements(EW_old))
EW_low = dblarr(n_elements(EW_old))
E_BV = dblarr(n_elements(lineSFR_old.E_BV))
age = dblarr(n_elements(age_old))
R_1 = [R_cutoff.xmm, R_cutoff.munics, R_cutoff.cosmos, R_cutoff.goods]
R = [R_cutoff.x_R, R_cutoff.m_R, R_cutoff.c_R, R_cutoff.g_R]

;remove, x, R_1

for i = 0, n_elements(R_1)-1 do begin
   if R_1[i] gt 0 then begin
      IDs[i] = IDs_old[i]
      mass[i] = mass_old[i]
      z[i] = z_old[i]
      lineSFR[i] = lineSFR_old.SFR[i]
      lineSFR_up[i] = lineSFR_old.SFR_up[i]
      lineSFR_low[i] = lineSFR_old.SFR_low[i]
      EW[i] = EW_old[i]
      EW_up[i] = EW_old_up[i]
      EW_low[i] = EW_old_low[i]
      E_BV[i] = lineSFR_old.E_BV[i]
      age[i] = age_old[i]
      lum[i] = lineSFR_old.lum[i]
      SFR_red[i] = lineSFR_old.SFR_red[i]
   endif
endfor


;index = where(z eq 0.0)
;remove, index, z, mass, lineSFR, lineSFR_up, lineSFR_low, IDs, EW,
;EW_up, EW_low, E_BV, age, E_BV_robin, lineSFR_robin, UVSFR_robin,
;Rmag, ID_robin, chi, lum

; IDs of X-ray sources that are AGN candidates - from Robin
AGN = [52, 73, 122, 239, 267, 366, 369, 410, 429, 464]
x = intarr(n_elements(AGN))
for i = 0, n_elements(AGN)-1 do begin
   x[i] = where(IDs eq AGN[i])
endfor

xray_mass = mass[x]
xray_lineSFR = lineSFR[x]
xray_IDs = IDs[x]               ; This line just gives back the AGN array, whatevs.
xray_lineSFR_up = lineSFR_up[x]
xray_lineSFR_low = lineSFR_low[x]
xray_z = z[x]
xray_EW = EW[x]
xray_EW_up = EW_up[x]
xray_EW_low = EW_low[x]
remove, x, mass, lineSFR, IDs, z, lineSFR_up, lineSFR_low, EW, EW_up, EW_low, E_BV, age, E_BV_robin, lineSFR_robin, UVSFR_robin, Rmag, ID_robin, chi, lum, SFR_red, xerr_up, xerr_low, R

; This section removes the galaxies with fewer than 5 photometric
; bands
four = [54, 70, 79, 81, 88, 97, 112]
x = intarr(n_elements(four))
for i = 0, n_elements(four)-1 do begin
   x[i] = where(IDs eq four[i])
endfor
remove, x, mass, lineSFR, IDs, z, lineSFR_up, lineSFR_low, EW, EW_up, EW_low, E_BV, age, E_BV_robin, lineSFR_robin, UVSFR_robin, Rmag, ID_robin, chi, lum, SFR_red, xerr_up, xerr_low, R

EW_old = EW
EW = EW/(1+z)
EW_up = EW_up/(1+z)
EW_low = EW_low/(1+z)
xray_EW = xray_EW/(1+xray_z)
xray_EW_up = xray_EW_up/(1+xray_z)
xray_EW_low = xray_EW_low/(1+xray_z)

;E_BV = garn_best(mass)

; This following section is for if not all galaxies made it
; through.  It basically pulls out the surviving z's from the
; cosmos_OII structure.  If all galaxies made it through, comment out
; this section
;ind = fix(strmid(goods_file, 30, 3))  ;Change '47' depending on length of file names!!!!!
;x = intarr(n_elements(ind))
;for i = 0, n_elements(x)-1 do begin
;   x[i] = where(all_OII.goods.id eq ind[i])
;endfor
;y = indgen(n_elements(all_OII.goods.id))
;remove, x, y;
;z_goods = all_OII.goods.z
;remove, y, z_goods
;z = [all_OII.xmm.z, all_OII.munics.z, all_OII.cosmos.z, z_goods]

; This section splits up the SFR and mass by redshift, for extra nice plotting
mass_lowz = dblarr(n_elements(mass))
mass_mid1z = dblarr(n_elements(mass))
mass_mid2z = dblarr(n_elements(mass))
mass_hiz =  dblarr(n_elements(mass))
lowz =  dblarr(n_elements(mass))
mid1z =  dblarr(n_elements(mass))
mid2z =  dblarr(n_elements(mass))
hiz =  dblarr(n_elements(mass))
lineSFR_lowz = dblarr(n_elements(mass))
lineSFR_mid1z = dblarr(n_elements(mass))
lineSFR_mid2z = dblarr(n_elements(mass))
lineSFR_hiz = dblarr(n_elements(mass))
lineSSFR_lowz = dblarr(n_elements(mass))
lineSSFR_mid1z = dblarr(n_elements(mass))
lineSSFR_mid2z = dblarr(n_elements(mass))
lineSSFR_hiz = dblarr(n_elements(mass))
EW_lowz = dblarr(n_elements(mass))
EW_mid1z = dblarr(n_elements(mass))
EW_mid2z = dblarr(n_elements(mass))
EW_hiz = dblarr(n_elements(mass))
E_BV_lowz = dblarr(n_elements(mass))
E_BV_mid1z = dblarr(n_elements(mass))
E_BV_mid2z = dblarr(n_elements(mass))
E_BV_hiz = dblarr(n_elements(mass))
lineSFR_robin_lowz = dblarr(n_elements(mass))
lineSFR_robin_mid1z = dblarr(n_elements(mass))
lineSFR_robin_mid2z = dblarr(n_elements(mass))
lineSFR_robin_hiz = dblarr(n_elements(mass))
UVSFR_robin_lowz = dblarr(n_elements(mass))
UVSFR_robin_mid1z = dblarr(n_elements(mass))
UVSFR_robin_mid2z = dblarr(n_elements(mass))
UVSFR_robin_hiz = dblarr(n_elements(mass))
lineSFR_red_lowz = dblarr(n_elements(mass))
lineSFR_red_mid1z = dblarr(n_elements(mass))
lineSFR_red_mid2z = dblarr(n_elements(mass))
lineSFR_red_hiz = dblarr(n_elements(mass))
lum_lowz = dblarr(n_elements(mass))
lum_mid1z= dblarr(n_elements(mass))
lum_mid2z = dblarr(n_elements(mass))
lum_hiz = dblarr(n_elements(mass))

for i = 0, n_elements(z)-1 do begin
   if z[i] ge 0 and z[i] lt 0.2 then begin
      mass_lowz[i] = mass[i]
      lineSFR_lowz[i] = lineSFR[i]
      lineSSFR_lowz[i] = lineSFR[i]-mass[i]
      lowz[i] = z[i]
      EW_lowz[i] = EW[i]
      E_BV_lowz[i] = E_BV[i]
      lineSFR_robin_lowz[i] = lineSFR_robin[i]
      UVSFR_robin_lowz[i] = UVSFR_robin[i]
      lineSFR_red_lowz[i] = SFR_red[i]
      lum_lowz[i] = lum[i]
   endif
   if z[i] ge 0.2 and z[i] lt 0.32 then begin
      mass_mid1z[i] = mass[i]
      lineSFR_mid1z[i] = lineSFR[i]
      lineSSFR_mid1z[i] = lineSFR[i]-mass[i]
      mid1z[i] = z[i]
      EW_mid1z[i] = EW[i]
      E_BV_mid1z[i] = E_BV[i]
      lineSFR_robin_mid1z[i] = lineSFR_robin[i]
      UVSFR_robin_mid1z[i] = UVSFR_robin[i]
      lineSFR_red_mid1z[i] = SFR_red[i]
      lum_mid1z[i] = lum[i]
   endif
   if z[i] ge 0.32 and z[i] lt 0.45 then begin
      mass_mid2z[i] = mass[i]
      lineSFR_mid2z[i] = lineSFR[i]
      lineSSFR_mid2z[i] = lineSFR[i]-mass[i]
      mid2z[i] = z[i]
      EW_mid2z[i] = EW[i]
      E_BV_mid2z[i] = E_BV[i]
      lineSFR_robin_mid2z[i] = lineSFR_robin[i]
      UVSFR_robin_mid2z[i] = UVSFR_robin[i]
      lineSFR_red_mid2z[i] = SFR_red[i]
      lum_mid2z[i] = lum[i]
   endif
   if z[i] ge 0.45 and z[i] le 0.57 then begin
      mass_hiz[i] = mass[i]
      lineSFR_hiz[i] = lineSFR[i]
      lineSSFR_hiz[i] = lineSFR[i]-mass[i]
      hiz[i] = z[i]
      EW_hiz[i] = EW[i]
      E_BV_hiz[i] = E_BV[i]
      lineSFR_robin_hiz[i] = lineSFR_robin[i]
      UVSFR_robin_hiz[i] = UVSFR_robin[i]
      lineSFR_red_hiz[i] = SFR_red[i]
      lum_hiz[i] = lum[i]
   endif
endfor

; Divide mass into chunks, find the median and +/- stddev for each bin
; for the purpose of plotting

; SFR calculation
bit1 = (lineSFR)[where(mass lt 8.15)]
bit2 = (lineSFR)[where(mass ge 8.15 and mass lt 8.75)]
bit3 = (lineSFR)[where(mass ge 8.75 and mass lt 9.35)]
bit4 = (lineSFR)[where(mass ge 9.35 and mass lt 9.95)]
bit5 = (lineSFR)[where(mass ge 9.95 and mass lt 10.55)]
bit6 = (lineSFR)[where(mass ge 10.55 and mass lt 11.05)]
bit7 = (lineSFR)[where(mass ge 11.05 and mass lt 11.85)]
;bit8 = (lineSFR)[where(mass ge 10.85 and mass lt 11.35)]
;bit9 = (lineSFR)[where(mass ge 11.35 and mass lt 11.85)]
;bit10 = (lineSFR)[where(mass ge 11.85 and mass lt 12.35)]

std1 = stddev(bit1)
std2 = stddev(bit2)
std3 = stddev(bit3)
std4 = stddev(bit4)
std5 = stddev(bit5)
std6 = stddev(bit6)
std7 = stddev(bit7)
;std8 = stddev(bit8)
;std9 = stddev(bit9)
;std10 = stddev(bit10)

med1 = median(bit1, /even)
med2 = median(bit2, /even)
med3 = median(bit3, /even)
med4 = median(bit4, /even)
med5 = median(bit5, /even)
med6 = median(bit6, /even)
med7 = median(bit7, /even)
;med8 = median(bit8, /even)
;med9 = median(bit9, /even)
;med10 = median(bit10, /even)


mean = [avg(bit1), avg(bit2), avg(bit3), avg(bit4), avg(bit5), avg(bit6), avg(bit7)]
med = [med1, med2, med3, med4, med5, med6, med7];, med8, med9, med10]
;std_up = [med1+std1, med2+std2, med3+std3, med4+std4, med5+std5, med6+std6, med7+std7];, med8+std8, med9+std9, med10+std10]
;std_low = [med1-std1, med2-std2, med3-std3, med4-std4, med5-std5, med6-std6, med7-std7];, med8-std8, med9-std9, med10-std10]
std_up = [mean[0]+std1, mean[1]+std2, mean[2]+std3, mean[3]+std4, mean[4]+std5, mean[5]+std6, mean[6]+std7]
std_low = [mean[0]-std1, mean[1]-std2, mean[2]-std3, mean[3]-std4, mean[4]-std5, mean[5]-std6, mean[6]-std7]
std = [std1, std2, std3, std4, std5, std6]
mass_ind = [7.85, 8.45, 9.05, 9.65, 10.25, 10.85, 11.45];, 11.1, 11.6, 12.1]


; SSFR calculation
other1 = (lineSFR-mass)[where(mass lt 8.15)]
other2 = (lineSFR-mass)[where(mass ge 8.15 and mass lt 8.85)]
other3 = (lineSFR-mass)[where(mass ge 8.85 and mass lt 9.55)]
other4 = (lineSFR-mass)[where(mass ge 9.55 and mass lt 10.25)]
other5 = (lineSFR-mass)[where(mass ge 10.25 and mass lt 10.95)]
other6 = (lineSFR-mass)[where(mass ge 10.95 and mass lt 11.8)]
;other7 = (lineSFR-mass)[where(mass ge 10.35 and mass lt 10.85)]
;other8 = (lineSFR-mass)[where(mass ge 10.85 and mass lt 11.35)]
;other9 = (lineSFR-mass)[where(mass ge 11.35 and mass lt 11.85)]
;other10 = (lineSFR-mass)[where(mass ge 11.85 and mass lt 12.35)]

ostd1 = stddev(other1)
ostd2 = stddev(other2)
ostd3 = stddev(other3)
ostd4 = stddev(other4)
ostd5 = stddev(other5)
ostd6 = stddev(other6)
;ostd7 = stddev(other7)
;ostd8 = stddev(other8)
;ostd9 = stddev(other9)
;ostd10 = stddev(other10)

omed1 = median(other1, /even)
omed2 = median(other2, /even)
omed3 = median(other3, /even)
omed4 = median(other4, /even)
omed5 = median(other5, /even)
omed6 = median(other6, /even)
;omed7 = median(other7, /even)
;omed8 = median(other8, /even)
;omed9 = median(other9, /even)
;omed10 = median(other10, /even)

omed = [omed1, omed2, omed3, omed4, omed5, omed6];, omed7, omed8, omed9];, omed10]
ostd_up = [omed1+ostd1, omed2+ostd2, omed3+ostd3, omed4+ostd4, omed5+ostd5, omed6+ostd6];, omed7+ostd7, omed8+ostd8, omed9+ostd9];, omed10+ostd10]
ostd_low = [omed1-ostd1, omed2-ostd2, omed3-ostd3, omed4-ostd4, omed5-ostd5, omed6-ostd6];, omed7-ostd7, omed8-ostd8, omed9-ostd9];, omed10-ostd10]

; Import Noeske data for comparison
readcol, 'data/PN_Mass_SFR/Kai0.2-0.45.txt', N_mass1, N_SFR1, N_mass_up1, N_SFR_up1, N_mass_low1, N_SFR_low1, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.45-0.70.txt', N_mass2, N_SFR2, N_mass_up2, N_SFR_up2, N_mass_low2, N_SFR_low2, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.70-0.85.txt', N_mass3, N_SFR3, N_mass_up3, N_SFR_up3, N_mass_low3, N_SFR_low3, format = 'D, D, D, D, D, D'
readcol, 'data/PN_Mass_SFR/Kai0.85-1.10.txt', N_mass4, N_SFR4, N_mass_up4, N_SFR_up4, N_mass_low4, N_SFR_low4, format = 'D, D, D, D, D, D'

; Import Pirzkal data for comparison
;P = pirzkal()
readcol, 'data/PN_Mass_SFR/SFR_Mass_bin/SFR_Mass_0.2_0.5.txt', P_mass1, P_SFR1, P_SFR_err1, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_bin/SFR_Mass_0.5_0.7.txt', P_mass2, P_SFR2, P_SFR_err2, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_bin/SFR_Mass_0.7_0.8.txt', P_mass3, P_SFR3, P_SFR_err3, format = 'D, D, D'
readcol, 'data/PN_Mass_SFR/SFR_Mass_bin/SFR_Mass_0.8_1.1.txt', P_mass4, P_SFR4, P_SFR_err4, format = 'D, D, D'

; liner fit to SFR v. M
fit = linfit(mass, lineSFR)
m = [4, 13]
s = fit[0]+m*fit[1]

;i = where(mass_lowz eq 0)
;remove, i, mass_lowz, lineSFR_lowz, lineSSFR_lowz, EW_lowz, E_BV_lowz, lineSFR_robin_lowz, UVSFR_robin_lowz, lineSFR_red_lowz
;i = where(mass_mid1z eq 0)
;remove, i, mass_mid1z, lineSFR_mid1z, lineSSFR_mid1z, EW_mid1z, E_BV_mid1z, lineSFR_robin_mid1z, UVSFR_robin_mid1z, lineSFR_red_mid1z
;i = where(mass_mid2z eq 0)
;remove, i, mass_mid2z, lineSFR_mid2z, lineSSFR_mid2z, EW_mid2z, E_BV_mid2z, lineSFR_robin_mid2z, UVSFR_robin_mid2z, lineSFR_red_mid2z
;i = where(mass_hiz eq 0)
;remove, i, mass_hiz, lineSFR_hiz, lineSSFR_hiz, EW_hiz, E_BV_hiz, lineSFR_robin_hiz, UVSFR_robin_hiz, lineSFR_red_hiz

; deredden the SFR 80% limit lines
wp = 0.3727
k_prime = 1.78 + (1.17 * ((0.011 / (wp^3)) - (0.198 / (wp^2)) + (1.509 / wp) - 2.156)) 
k_prime = k_prime/0.44    ; As per Calzetti for emission line extinction
E_stars_low = 0.21404000      ;median(E_BV_lowz, /even) ;0.3030657
E_stars_mid1 = 0.23430880    ;median(E_BV_mid1z, /even) ;0.38563870
E_stars_mid2 = 0.26536277    ;median(E_BV_mid2z, /even) ;0.25885230
E_stars_hi = 0.27781999         ;median(E_BV_hiz, /even) ;0.30247930
log_lum_low = 10^(0.4*E_stars_low*k_prime)
log_lum_mid1 = 10^(0.4*E_stars_mid1*k_prime)
log_lum_mid2 = 10^(0.4*E_stars_mid2*k_prime)
log_lum_hi = 10^(0.4*E_stars_hi*k_prime)
E_stars_mid = 0.26085001      ;median([E_BV_mid1z, E_BV_mid2z], /even)
log_lum_mid = 10^(0.4*E_stars_mid*k_prime)
E_stars = median(E_BV, /even)
log_lum = 10^(0.4*E_stars*k_prime)


;print, 41.0509+alog10(log_lum_hi) - 41.1817d
;print, 40.7627+alog10(log_lum_mid2) - 41.1817d
;print, 40.4329+alog10(log_lum_mid1) - 41.1817d
;print,  39.6548+alog10(log_lum_low) - 41.1817d


lowz80 = 39.5726+alog10(log_lum_low) - 41.1817d
mid1z80 = 40.3381+alog10(log_lum_mid1) - 41.1817d
mid2z80 = 40.6794+alog10(log_lum_mid2) - 41.1817d
hiz80 = 40.958+alog10(log_lum_hi) - 41.1817d
new_mass = [mass_lowz[where(mass_lowz gt 7)], mass_mid1z[where(mass_mid1z gt 9.5)], mass_mid2z[where(mass_mid2z gt 10)], mass_hiz[where(mass_hiz gt 10.5)]]
new_lineSFR = [lineSFR_lowz[where(mass_lowz gt 7)], lineSFR_mid1z[where(mass_mid1z gt 9.5)], lineSFR_mid2z[where(mass_mid2z gt 10)], lineSFR_hiz[where(mass_hiz gt 10.5)]] 
for k = 0, 100 do begin
   fit = linfit(new_mass, new_lineSFR, sigma = sig)
   m = fit[1]
   b = fit[0]
   scatter = stdev(new_lineSFR - (m*new_mass + b))
   newline = mass*m+b
   cut1 = (lowz80-b)/m
   cut2 = (mid1z80-b)/m
   cut3 = (mid2z80-b)/m
   cut4 = (hiz80-b)/m
   ;print, cut1, cut2, cut3, cut4
   new_mass = [mass_lowz[where(mass_lowz gt cut1)], mass_mid1z[where(mass_mid1z gt cut2)], mass_mid2z[where(mass_mid2z gt cut3)], mass_hiz[where(mass_hiz gt cut4)]]
   new_lineSFR = [lineSFR_lowz[where(mass_lowz gt cut1)], lineSFR_mid1z[where(mass_mid1z gt cut2)], lineSFR_mid2z[where(mass_mid2z gt cut3)], lineSFR_hiz[where(mass_hiz gt cut4)]]
   k += 1
endfor

fit = linfit(new_mass, new_lineSFR, sigma = sig)

;print, mass_lowz[where(mass_lowz gt cut1, complement = g)], lineSFR_lowz[where(mass_lowz gt cut1)]

; PLOT ALL THE THINGS!!! Note that I have to remember to have moved the
; old plots into a different folder, or they will be written over with
; these
;A = FINDGEN(24) *  (!PI*2/24.) 
;USERSYM, COS(A), SIN(A), /FILL 
loadct, 17
openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSFRvsmass_z'
plot, mass_lowz, lineSFR, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7.45, 12], ystyle = 1, yrange = [-2, 3.5]
hline, 40.958+alog10(log_lum_hi) - 41.1817d, lines=1, color = 55 ; 80% completeness at z = 0.51
hline, 40.6794+alog10(log_lum_mid2) - 41.1817d, lines=1, color = 140 ; 80% completeness at z = 0.4
hline, 40.3381+alog10(log_lum_mid1) - 41.1817d, lines=1, color = 100 ; 80% completeness at z = 0.275
hline, 39.5726+alog10(log_lum_low) - 41.1817d, lines=1, color = 242 ; 80% completeness at z = 0.1
;oplot, mass_lowz, lineSFR, psym = sym(1), symsize = 1.2, color =  242
;oplot, mass_mid1z,lineSFR, psym = sym(1), symsize = 1.2, color = 100
;oplot, mass_mid2z, lineSFR, psym = sym(1), symsize = 1.2, color = 140
;oplot, mass_hiz, lineSFR, psym = sym(1), symsize = 1.2, color = 55

oplot, mass_lowz[where(mass_lowz gt cut1, complement = g)], lineSFR_lowz[where(mass_lowz gt cut1)], color = 242, symsize = 1.1, psym = sym(1)
oplot, mass_lowz[g], lineSFR_lowz[g], psym = sym(6), color = 242, symsize = 1.1, thick = 1.5     ;2

oplot, mass_mid1z[where(mass_mid1z gt cut2, complement = g)], lineSFR_mid1z[where(mass_mid1z gt cut2)], color = 100, symsize = 1.1, psym = sym(1)
oplot, mass_mid1z[g], lineSFR_mid1z[g], psym =sym(6), color = 100, symsize = 1.1  ;6

oplot, mass_mid2z[where(mass_mid2z gt cut3, complement = g)], lineSFR_mid2z[where(mass_mid2z gt cut3)], color = 140, symsize = 1.1, psym = sym(1)
oplot, mass_mid2z[g], lineSFR_mid2z[g], psym =sym(6), color = 140, symsize = 1.1  ;5

oplot, mass_hiz[where(mass_hiz gt cut4, complement = g)], lineSFR_hiz[where(mass_hiz gt cut4)], color = 55, symsize = 1.1, psym = sym(1)
oplot, mass_hiz[g], lineSFR_hiz[g], psym = sym(6), color = 55, symsize = 1.1             ;4



;oplot, mass_ind, mean, psym = sym(1), symsize = 1.9
;oplot, mass_ind, std_up, lines=3
;oplot, mass_ind, std_low, lines=3
legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [sym(1), sym(1), sym(1), sym(1)], color = [242, 100, 140, 45]
oploterror, [11.6], [-1.55], [alog10(median(10^xerr_low))], [0], /lobar      
oploterror, [11.6], [-1.55], [alog10(median(10^xerr_up))], [0], /hibar 
oplot, [11.6], [-1.55], psym = 3
oploterror, [11.6], [-1.55], [0], [alog10(median(10^lineSFR_up))], /hibar
oploterror, [11.6], [-1.55], [0], [alog10(median(10^lineSFR_low))], /lobar
oplot, xray_mass[[1, 3, 5]], xray_lineSFR[[1, 3, 5]], psym = 7, symsize = 2, color = 242
oplot, xray_mass[[2, 4, 8]], xray_lineSFR[[2, 4, 8]], psym = 7, symsize = 2, color = 100
oplot, [xray_mass[0]], [xray_lineSFR[0]], psym = 7, symsize = 2, color = 140
oplot, xray_mass[[6, 7, 9]], xray_lineSFR[[6, 7, 9]], psym = 7, symsize = 2, color = 55
;oplot, m, s, thick = 4 ; Linear fit to the data, I don't think I'll include it
closepps


new_lum = [lum_lowz[where(mass_lowz gt cut1)], lum_mid1z[where(mass_mid1z gt cut2)], lum_mid2z[where(mass_mid2z gt cut3)], lum_hiz[where(mass_hiz gt cut4)]] 

loadct, 17
openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSFRvsmass_masscut'
plot, new_mass, new_lineSFR, psym =3,  xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7, 12], ystyle = 1, yrange = [-2, 2.5]
oplot, new_mass, new_lineSFR, psym = sym(1), symsize = 1.2, color = 242
oploterror, [11.6], [-1.55], [alog10(median(10^xerr_low))], [0], /lobar      
oploterror, [11.6], [-1.55], [alog10(median(10^xerr_up))], [0], /hibar 
oploterror, [11.6], [-1.55], [0], [alog10(median(10^lineSFR_up))], /hibar
oploterror, [11.6], [-1.55], [0], [alog10(median(10^lineSFR_low))], /lobar
oplot, [7.4, 11.7], fit[0]+fit[1]*[7.4, 11.7], lines = 4
closepps


; Only do this for the mass-cut sample
met = -10.8297+3.6478*new_mass-0.16706*new_mass^2   ; Lara-Lopez mass metallicity relationship
new_SFR = alog10(7.9d-42*new_lum/(-1.75*met + 16.73))     ; new SFR calculation using metallicity

bit1 = (new_lineSFR)[where(new_mass lt 8.15)]
bit2 = (new_lineSFR)[where(new_mass ge 8.15 and new_mass lt 8.75)]
bit3 = (new_lineSFR)[where(new_mass ge 8.75 and new_mass lt 9.35)]
bit4 = (new_lineSFR)[where(new_mass ge 9.35 and new_mass lt 9.95)]
bit5 = (new_lineSFR)[where(new_mass ge 9.95 and new_mass lt 10.55)]
bit6 = (new_lineSFR)[where(new_mass ge 10.55 and new_mass lt 11.05)]
bit7 = (new_lineSFR)[where(new_mass ge 11.05 and new_mass lt 11.85)]

med1 = median(bit1, /even)
med2 = median(bit2, /even)
med3 = median(bit3, /even)
med4 = median(bit4, /even)
med5 = median(bit5, /even)
med6 = median(bit6, /even)
med7 = median(bit7, /even)

med = [med1, med2, med3, med4, med5, med6, med7]

mbit1 = (new_SFR)[where(new_mass lt 8.15)]
mbit2 = (new_SFR)[where(new_mass ge 8.15 and new_mass lt 8.75)]
mbit3 = (new_SFR)[where(new_mass ge 8.75 and new_mass lt 9.35)]
mbit4 = (new_SFR)[where(new_mass ge 9.35 and new_mass lt 9.95)]
mbit5 = (new_SFR)[where(new_mass ge 9.95 and new_mass lt 10.55)]
mbit6 = (new_SFR)[where(new_mass ge 10.55 and new_mass lt 11.05)]
mbit7 = (new_SFR)[where(new_mass ge 11.05 and new_mass lt 11.85)]

mstd1 = stddev(mbit1)
mstd2 = stddev(mbit2)
mstd3 = stddev(mbit3)
mstd4 = stddev(mbit4)
mstd5 = stddev(mbit5)
mstd6 = stddev(mbit6)
mstd7 = stddev(mbit7)

mmed1 = median(mbit1, /even)
mmed2 = median(mbit2, /even)
mmed3 = median(mbit3, /even)
mmed4 = median(mbit4, /even)
mmed5 = median(mbit5, /even)
mmed6 = median(mbit6, /even)
mmed7 = median(mbit7, /even)

mmed = [mmed1, mmed2, mmed3, mmed4, mmed5, mmed6, mmed7]
mstd = [mstd1, mstd2, mstd3, mstd4, mstd5, mstd6, mstd7]
mstd_up = [mmed1+mstd1, mmed2+mstd2, mmed3+mstd3, mmed4+mstd4, mmed5+mstd5, mmed6+mstd6, mmed7+mstd7]
mstd_low = [mmed1-mstd1, mmed2-mstd2, mmed3-mstd3, mmed4-mstd4, mmed5-mstd5, mmed6-mstd6, mmed7-mstd7]

;Do some bootstrapping
b1 = stdev(bbootstrap(bit1, 200, funct = 'median'))
b2 = stdev(bbootstrap(bit2, 200, funct = 'median'))
b3 = stdev(bbootstrap(bit3, 200, funct = 'median'))
b4 = stdev(bbootstrap(bit4, 200, funct = 'median'))
b5 = stdev(bbootstrap(bit5, 200, funct = 'median'))
b6 = stdev(bbootstrap(bit6, 200, funct = 'median'))
b7 = stdev(bbootstrap(bit7, 200, funct = 'median'))

b = [b1, b2, b3, b4, b5, b6, b7]

mb1 = stdev(bbootstrap(mbit1, 200, funct = 'median'))
mb2 = stdev(bbootstrap(mbit2, 200, funct = 'median'))
mb3 = stdev(bbootstrap(mbit3, 200, funct = 'median'))
mb4 = stdev(bbootstrap(mbit4, 200, funct = 'median'))
mb5 = stdev(bbootstrap(mbit5, 200, funct = 'median'))
mb6 = stdev(bbootstrap(mbit6, 200, funct = 'median'))
mb7 = stdev(bbootstrap(mbit7, 200, funct = 'median'))

mb = [mb1, mb2, mb3, mb4, mb5, mb6, mb7]

line1 = linfit(new_mass, new_lineSFR, sigma = sig1)
line2 = linfit(new_mass, new_SFR, sigma = sig2)

; Redo the above plot but with new SFR
openpps, 'SFR_plots/HPS_o2sed_IRAC/metallicity_SFR'
plot, mass_ind, med, psym = sym(1), symsize = 1.9, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7.5, 12], ystyle = 1, yrange = [-1.5, 1.75]
oplot, [7,13], line1[0]+line1[1]*[7,13], lines = 3
oploterror, mass_ind, med, b;, /hibar 
oplot, [7,13], line2[0]+line2[1]*[7,13], color = 140, lines = 3
oplot, mass_ind, mmed, psym = sym(1), symsize = 1.9, color = 140
oploterror, mass_ind, mmed, mb, color = 140;, /lobar   
closepps

print, 'Old line fit:', line1[0], line1[1], sig1
print, 'New metallicity line fit:', line2[0], line2[1], sig2

; Lara-Lopez z out to 0.1 local galaxy line
m = [4, 13]
s = -5.3126+0.5547*m

openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSFRvsmass_z_PN'
!p.multi = [0, 3, 2]
plot, mass_lowz, lineSFR, psym = sym(1), symsize = .7, ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [6, 12], ystyle = 1, yrange = [-2, 3.7], pos = [.05, .5, .35, 0.95], XTickformat='(A1)', ytickname = [' ', '-1', '0', '1', '2', '3'], charsize = 1.75
oplot, m, s, color = 242, thick = 4
oplot, m, s+0.349, color = 242, thick = 4, lines = 3
oplot, m, s-0.349, color = 242, thick = 4, lines = 3
oplot, mass_lowz, lineSFR, psym = sym(1), symsize = .8, color = 55
oplot, xray_mass[[1, 3, 5]], xray_lineSFR[[1, 3, 5]], psym = 7, symsize = 1.5, color = 55
xyouts, 9.5, 3.2,  '0.0 < z < 0.2', charsize = .9
hline, 39.5726+alog10(log_lum_low) - 41.1817d, lines=2, color = 55 ; 80% completeness at z = 0.1
legend, ['This survey'], psym = [sym(1)], charsize = .75, color=[55]


plot, mass_mid1z, lineSFR, psym = sym(1), symsize = .7, xstyle = 1, xrange = [6, 12], ystyle = 1, yrange = [-2, 3.7], pos = [.35, .5, .65, .95], XTickformat='(A1)', YTickformat='(A1)'
oplot, m, s, color = 242, thick = 4
oplot, m, s+0.349, color = 242, thick = 4, lines = 3
oplot, m, s-0.349, color = 242, thick = 4, lines = 3
oplot, mass_mid1z, lineSFR, psym = sym(1), symsize = .8, color = 55
oplot, xray_mass[[2, 4, 8]], xray_lineSFR[[2, 4, 8]], psym = 7, symsize = 1.5, color = 55
oplot, [xray_mass[0]], [xray_lineSFR[0]], psym = 7, symsize = 1.5, color = 55
hline, 40.4823+alog10(log_lum_mid) - 41.1817d, lines=2, color = 55 ; 80% completeness at z = 0.325
oplot, mass_mid2z, lineSFR, psym = sym(1), symsize = .8, color = 55
oplot, N_mass1, N_SFR1, psym = sym(1), color = 100
oplot, N_mass1, N_SFR_low1, lines = 3, color = 100
oplot, N_mass1, N_SFR_up1, lines=3, color = 100
oplot, P_mass1, P_SFR1, psym = sym(1), color = 140
oplot, P_mass1, P_SFR1+P_SFR_err1, lines=3, color=140
oplot, p_mass1, P_SFR1-P_SFR_err1, lines=3, color=140
xyouts, 9.5, 3.2,  '0.2 < z < 0.45', charsize = .9
legend, ['This survey',  'Noeske et al. 2007', 'Pirzkal et al. 2012'], psym = [sym(1), sym(1), sym(1)], color = [55, 100, 140], charsize = .75


plot, mass_hiz, lineSFR, psym = sym(1), symsize = .7, xtitle = textoidl('log(M_* [M_{sun}])'), xstyle = 1, xrange = [6, 12], ystyle = 1, yrange = [-2, 3.7], pos = [.65, .5, .95, .95], YTickformat='(A1)', xtickname=[' ', '7', '8', '9', '10', '11', '12'], charsize = 1.75
oplot, m, s, color = 242, thick = 4
oplot, m, s+0.349, color = 242, thick = 4, lines = 3
oplot, m, s-0.349, color = 242, thick = 4, lines = 3
oplot, mass_hiz, lineSFR, psym = sym(1), symsize = .8, color = 55
oplot, xray_mass[[6, 7, 9]], xray_lineSFR[[6, 7, 9]], psym = 7, symsize = 1.5, color = 55
hline, 40.9516+alog10(log_lum_hi) - 41.1817d, lines=2, color = 55               ; 80% completeness at z = 0.505
oplot, N_mass2, N_SFR2, psym = sym(1), color = 100
oplot, N_mass2, N_SFR_low2, lines = 3, color = 100
oplot, N_mass2, N_SFR_up2, lines=3, color = 100
oplot, P_mass2, P_SFR2, psym = sym(1), color = 140
oplot, P_mass2, P_SFR2+P_SFR_err2, lines=3, color=140
oplot, p_mass2, P_SFR2-P_SFR_err2, lines=3, color=140
xyouts, 9.5, 3.2,  '0.45 < z < 0.7', charsize = .9
legend, [ 'This survey', 'Noeske et al. 2007', 'Pirzkal et al. 2012'], psym = [ sym(1), sym(1), sym(1)], color = [55, 100, 140], charsize = .75


plot, N_mass3, N_SFR3, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [6, 12], ystyle = 1, yrange = [-2, 3.7], pos = [.05, .1, .35, .5], xtickname=['6', '7', '8', '9', '10', '11', ' '], ytickname = ['-2', '-1', '0', '1', '2', '3'], charsize = 1.75
oplot, m, s, color = 242, thick = 4
oplot, m, s+0.349, color = 242, thick = 4, lines = 3
oplot, m, s-0.349, color = 242, thick = 4, lines = 3
oplot, N_mass3, N_SFR3, psym = sym(1), color = 100
xyouts, 9.5, 3.2,  '0.7 < z < 0.85', charsize = .9
oplot, N_mass3, N_SFR_low3, lines = 3, color = 100
oplot, N_mass3, N_SFR_up3, lines=3, color = 100
oplot, P_mass3, P_SFR3, psym = sym(1), color = 140
oplot, P_mass3, P_SFR3+P_SFR_err3, lines=3, color=140
oplot, p_mass3, P_SFR3-P_SFR_err3, lines=3, color=140
legend, ['Noeske et al. 2007', 'Pirzkal et al. 2012'], psym = [ sym(1), sym(1)], color = [ 100, 140], charsize = .75


plot, N_mass4, N_SFR4, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), xstyle = 1, xrange = [6, 12], ystyle = 1, yrange = [-2, 3.7], pos = [.35, .1, .65, .5] , YTickformat='(A1)', xtickname=[' ', '7', '8', '9', '10', '11', '12'], charsize = 1.75
oplot, m, s, color = 242, thick = 4
oplot, m, s+0.349, color = 242, thick = 4, lines = 3
oplot, m, s-0.349, color = 242, thick = 4, lines = 3
oplot, N_mass4, N_SFR4, color = 100, psym = sym(1)
xyouts, 9.5, 3.2,  '0.85 < z < 1.1', charsize = .9
oplot, N_mass4, N_SFR_low4, lines = 3, color = 100
oplot, N_mass4, N_SFR_up4, lines=3, color = 100
oplot, P_mass4, P_SFR4, psym = sym(1), color = 140
oplot, P_mass4, P_SFR4+P_SFR_err4, lines=3, color=140
oplot, p_mass4, P_SFR4-P_SFR_err4, lines=3, color=140
legend, ['Noeske et al. 2007', 'Pirzkal et al. 2012'], psym = [sym(1), sym(1)], color = [ 100, 140], charsize = .75

closepps

; For to plot lines of constant SF
x = [10^6., 10^7., 10^8., 10^9., 10^10., 10^11., 10^12., 10^13., 10^14.]

;Propagate error to find SSFR error
;Had to do some fancy finagling because of the logginess
sig_ssfr_up = sqrt((10^(lineSFR-mass))^2*(((10^lineSFR_up)/(10^lineSFR))^2 + ((10^xerr_up)/(10^mass))^2))
sig_ssfr_low = sqrt((((10^lineSFR)/(10^mass))^2)*(((10^lineSFR_low)/(10^lineSFR))^2 + ((10^xerr_low)/(10^mass))^2))
ssfr_up = alog10(10^(-10.4d)*(median(sig_ssfr_up)/10^median(lineSFR-mass))+10^(-10.4d))+10.4d
ssfr_low = alog10(10^(-10.4d)*(median(sig_ssfr_low)/10^median(lineSFR-mass))+10^(-10.4d))+10.4d

openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSSFRvsmass_z'
plot, mass_lowz, lineSSFR_lowz, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(sSFR [yr^{-1}])'), xstyle = 1, xrange = [7.2, 12], ystyle = 1, yrange =[-11., -7]
oplot, alog10(x), -2-alog10(x), lines = 1
oplot, alog10(x), -alog10(x), lines = 1
oplot, alog10(x), 2-alog10(x), lines = 1
xyouts, 7.3, -9.2, textoidl('0.01 M_{sun} yr^{-1}'), charsize = 1.1
xyouts, 7.3, -7.2, textoidl('1 M_{sun} yr^{-1}'), charsize = 1.1
xyouts, 9.3, -7.2, textoidl('100 M_{sun} yr^{-1}'), charsize = 1.1
;hline, 40.958 - 41.1817d - median(mass_hiz_temp[sort(mass_hiz_temp)], /even), lines=1, color = 45 ; 80% completeness at z = 0.51
;hline, 40.6794 - 41.1817d - median(mass_mid2z_temp[sort(mass_mid2z_temp)],/even),  lines=1, color = 140 ; 80% completeness at z = 0.4
;hline, 40.3381 - 41.1817d - median(mass_mid1z_temp[sort(mass_mid1z_temp)], /even), lines=1, color = 100 ; 80% completeness at z = 0.275
;hline, 39.5726 - 41.1817d -median(mass_lowz_temp[sort(mass_lowz_temp)], /even), lines=1, color =242 ; 80% completeness at z = 0.1
oplot, mass_lowz[where(mass_lowz gt cut1, complement = g)], lineSSFR_lowz[where(mass_lowz gt cut1)], color = 242, symsize = 1.1, psym = sym(1)
oplot, mass_lowz[g], lineSSFR_lowz[g], psym = sym(6), color = 242, symsize = 1.1, thick = 1.5     ;2

oplot, mass_mid1z[where(mass_mid1z gt cut2, complement = g)], lineSSFR_mid1z[where(mass_mid1z gt cut2)], color = 100, symsize = 1.1, psym = sym(1)
oplot, mass_mid1z[g], lineSSFR_mid1z[g], psym =sym(6), color = 100, symsize = 1.1  ;6

oplot, mass_mid2z[where(mass_mid2z gt cut3, complement = g)], lineSSFR_mid2z[where(mass_mid2z gt cut3)], color = 140, symsize = 1.1, psym = sym(1)
oplot, mass_mid2z[g], lineSSFR_mid2z[g], psym =sym(6), color = 140, symsize = 1.1  ;5

oplot, mass_hiz[where(mass_hiz gt cut4, complement = g)], lineSSFR_hiz[where(mass_hiz gt cut4)], color = 55, symsize = 1.1, psym = sym(1)
oplot, mass_hiz[g], lineSSFR_hiz[g], psym = sym(6), color = 55, symsize = 1.1             ;4

;oplot, mass_ind, omed, psym = sym(1), symsize = 1.9
;oplot, mass_ind, ostd_up, lines=3
;oplot, mass_ind, ostd_low, lines=3
legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [sym(1), sym(1), sym(1), sym(1)], color = [242, 100, 140, 55], pos = [10.8, -7.12]
oploterror, [7.6], [-10.4], [alog10(median(10^xerr_up))], [0], /lobar      
oploterror, [7.6], [-10.4], [alog10(median(10^xerr_low))], [0], /hibar 
oplot, [7.6], [-0.4], psym = 3
oploterror, [7.6], [-10.4], [0], [(ssfr_up)], /hibar
oploterror, [7.6], [-10.4], [0], [(ssfr_low)], /lobar
oplot, xray_mass[[1, 3, 5]], xray_lineSFR[[1, 3, 5]]- xray_mass[[1, 3, 5]], psym = 7, symsize = 2, color = 242
oplot, xray_mass[[2, 4, 8]], xray_lineSFR[[2, 4, 8]]-xray_mass[[2, 4, 8]], psym = 7, symsize = 2, color = 100
oplot, [xray_mass[0]], [xray_lineSFR[0]-xray_mass[0]], psym = 7, symsize = 2, color = 140
oplot, xray_mass[[6, 7, 9]], xray_lineSFR[[6, 7, 9]]-xray_mass[[6, 7, 9]], psym = 7, symsize = 2, color = 55
closepps

;openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSSFRvsmass_masscut'
;plot, new_mass, new_lineSFR-new_mass, psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(sSFR [yr^{-1}])'), xstyle = 1, xrange = [7.2, 12], ystyle = 1, yrange =[-11., -7]
;oplot, new_mass, new_lineSFR-new_mass, psym = sym(1), color = 242, symsize = 1.1      ;2
;oplot, alog10(x), -2-alog10(x), lines = 1
;oplot, alog10(x), -alog10(x), lines = 1
;oplot, alog10(x), 2-alog10(x), lines = 1
;xyouts, 7.3, -9.2, textoidl('0.01 M_{sun} yr^{-1}'), charsize = 1.1
;xyouts, 7.3, -7.2, textoidl('1 M_{sun} yr^{-1}'), charsize = 1.1
;xyouts, 9.3, -7.2, textoidl('100 M_{sun} yr^{-1}'), charsize = 1.1
;oplot, mass_ind, omed, psym = sym(1), symsize = 1.9
;oplot, mass_ind, ostd_up, lines=3
;oplot, mass_ind, ostd_low, lines=3
;oploterror, [7.6], [-10.4], [alog10(median(10^xerr_up))], [0], /lobar      
;oploterror, [7.6], [-10.4], [alog10(median(10^xerr_low))], [0], /hibar 
;oplot, [7.6], [-0.4], psym = 3
;oploterror, [7.6], [-10.4], [0], [(ssfr_up)], /hibar
;oploterror, [7.6], [-10.4], [0], [(ssfr_low)], /lobar
;closepps

print, 'mass error up:',  alog10(median(10^xerr_up))  
print, 'mass error low:', alog10(median(10^xerr_low))
print, ssfr_up
print, ssfr_low
print, alog10(median(10^lineSFR_up))
print, alog10(median(10^lineSFR_low))

i = where(mass_lowz eq 0)
remove, i, mass_lowz, lineSFR_lowz, lineSSFR_lowz, EW_lowz, E_BV_lowz, lineSFR_robin_lowz, UVSFR_robin_lowz, lineSFR_red_lowz
i = where(mass_mid1z eq 0)
remove, i, mass_mid1z, lineSFR_mid1z, lineSSFR_mid1z, EW_mid1z, E_BV_mid1z, lineSFR_robin_mid1z, UVSFR_robin_mid1z, lineSFR_red_mid1z
i = where(mass_mid2z eq 0)
remove, i, mass_mid2z, lineSFR_mid2z, lineSSFR_mid2z, EW_mid2z, E_BV_mid2z, lineSFR_robin_mid2z, UVSFR_robin_mid2z, lineSFR_red_mid2z
i = where(mass_hiz eq 0)
remove, i, mass_hiz, lineSFR_hiz, lineSSFR_hiz, EW_hiz, E_BV_hiz, lineSFR_robin_hiz, UVSFR_robin_hiz, lineSFR_red_hiz



openpps, 'SFR_plots/HPS_o2sed_IRAC/mass_hist'
multiplot, [4, 1], /square
cgHistoplot, [mass_lowz, mass_mid1z], bin = 0.25, /line_fill, color = 'firebrick', thick = 6, polycolor = 'firebrick', orientation = [45, -45], xrange = [7,12.5],  yrange = [0, 25], xtickname=[  '7', '8', '9', '10', '11', '12', ' '], xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black', ytitle = 'Number of [O II] Emitters'
cgText, 7.3, 21.6, '0.0 < z 0.35', charsize = .9
multiplot
cgHistoplot, mass_mid2z,  bin = 0.25, /line_fill, color = 'firebrick', thick = 6, polycolor = 'firebrick', orientation = [45, -45], xrange = [7,12.5],  yrange = [0, 25], YTickformat='(A1)', xtickname=[ ' ', '8', '9', '10', '11', '12', ' '],  ytitle = ' ', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black' , xtitle = textoidl('log(M_* [M_{sun}])')
cgText, 7.3, 21.6, '0.35 < z 0.45', charsize = .9
multiplot
cgHistoplot, mass_hiz,  bin = 0.25, /line_fill, color = 'firebrick', thick = 6, polycolor = 'firebrick', orientation = [45, -45], xrange = [7,12.5],  yrange = [0, 25],YTickformat='(A1)',  xtickname=[ ' ', '8', '9', '10', '11', '12'], ytitle = ' ',xcharsize = '.75', ycharsize = '.75',  /oprobability, probcolor = 'black'
cgText, 7.3, 21.6, '0.45 < z 0.57', charsize = .9
axis, yaxis = 1, ycharsize = '.8', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks=5, ytitle = 'Cumulative Probability'
multiplot
closepps


openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSRF_hist'
multiplot, [4, 1], /square
cgHistoplot, [lineSFR_lowz, lineSFR_mid1z], bin = 0.2, /line_fill, color = 'navy', thick = 6, polycolor = 'navy', orientation = [45, -45], xrange = [-1.5, 2.5],  yrange = [0, 25],  ytitle = 'Number of [OII] Emitters', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black' , xtickname=[ '-1', '0', '1', '2']
cgText, -1.3, 21.3, '0.0 < z 0.35', charsize = .9
multiplot
cgHistoplot, lineSFR_mid2z,  bin = 0.2, /line_fill, color = 'navy', thick = 6, polycolor ='navy', orientation = [45, -45], xrange = [-1.5, 2.5],  yrange = [0, 25], YTickformat='(A1)',  ytitle = ' ', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black', xtickname=[ '-1', '0', '1', '2'],  xtitle = textoidl('log(SFR [M_{sun} yr^{-1}])')
cgText, -1.3, 21.3, '0.35 < z 0.45', charsize = .9
multiplot
cgHistoplot, lineSFR_hiz,  bin = 0.2, /line_fill, color = 'navy', thick = 6, polycolor = 'navy', orientation = [45, -45], xrange = [-1.5, 2.5],  yrange = [0, 25],YTickformat='(A1)',  ytitle = ' ',xcharsize = '.75', ycharsize = '.75',  /oprobability, probcolor = 'black', xtickname=[ '-1', '0', '1', '2']
cgText, -1.3, 21.3, '0.45 < z 0.57', charsize = .9
axis, yaxis = 1, ytitle = 'Cumulative Probability', ycharsize = '.8', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks = 5
closepps

;openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSRF_red_hist'
;multiplot, [4, 1], /square
;cgHistoplot, [lineSFR_red_lowz, lineSFR_red_mid1z], bin = 0.35, /line_fill, color = 'navy', thick = 6, polycolor = 'navy', orientation = [45, -45], xrange = [-1.5, 3.5],  yrange = [0, 40],  ytitle = 'Number of [OII] Emitters', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black' , xtickname=[ '-2', '-1', '0', '1', '2',  '3']
;cgText, 1.5, 22.6, '0.0 < z 0.35', charsize = .9
;multiplot
;cgHistoplot, lineSFR_red_mid2z,  bin = 0.35, /line_fill, color = 'navy', thick = 6, polycolor ='navy', orientation = [45, -45], xrange = [-2, 3.5],  yrange = [0, 40], YTickformat='(A1)',  ytitle = ' ', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black', xtickname=[ ' ', '-1', '0', '1', '2',  '3'],  xtitle = textoidl('log(SFR [M_{sun} yr^{-1}])')
;cgText, -1.7, 22.6, '0.35 < z 0.45', charsize = .9
;multiplot
;cgHistoplot, lineSFR_red_hiz,  bin = 0.35, /line_fill, color = 'navy', thick = 6, polycolor = 'navy', orientation = [45, -45], xrange = [-2, 3.5],  yrange = [0, 40],YTickformat='(A1)',  ytitle = ' ',xcharsize = '.75', ycharsize = '.75',  /oprobability, probcolor = 'black', xtickname=[ ' ', '-1', '0', '1', '2',  '3']
;cgText, -1.7, 22.6, '0.45 < z 0.57', charsize = .9
;axis, yaxis = 1, ytitle = 'Cumulative Probability', ycharsize = '.8', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks = 5
;closepps

openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSSRF_hist'
multiplot, [4, 1], /square, mxtitoffset = 1.1
cgHistoplot, [lineSSFR_lowz, lineSSFR_mid1z], bin = 0.25, /line_fill, color = 'forest green', thick = 6, polycolor = 'forest green', orientation = [45, -45], xrange = [-11, -6.5],  yrange = [0, 20], xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black' , xtickname=[ '-11', '-10', '-9', '-8',  '-7'],  ytitle = 'Number of [OII] Emitters'
cgText, -10.7, 17.7, '0.0 < z 0.35', charsize = .9
multiplot
cgHistoplot, lineSSFR_mid2z,  bin = 0.25, /line_fill, color =  'forest green', thick = 6, polycolor = 'forest green', orientation = [45, -45], xrange = [-11, -6.5],  yrange = [0, 20], YTickformat='(A1)', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black', xtickname=[ ' ', '-10',  '-9', '-8', '-7'], xtitle = textoidl('log(sSFR [yr^{-1}])'), ytitle = ' '
cgText, -10.7, 17.7, '0.35 < z 0.45', charsize = .9
multiplot
cgHistoplot, lineSSFR_hiz,  bin = 0.25, /line_fill, color = 'forest green', thick = 6, polycolor =  'forest green', orientation = [45, -45], xrange = [-11, -6.5],  yrange = [0, 20],YTickformat='(A1)',  ytitle = ' ',xcharsize = '.75', ycharsize = '.75',  /oprobability, probcolor = 'black', xtickname=[ ' ', '-10', '-9', '-8',  '-7']
cgText, -10.7, 17.7, '0.45 < z 0.57', charsize = .9
axis, yaxis = 1, ytitle = 'Cumulative Probability', ycharsize = '.8', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks = 5
closepps


;stop


readcol, 'data/robin_completeness.dat', complete_z, complete_90, complete_80, complete_70, complete_60, complete_50, format = 'D,D,D,D,D,D'
complete_90_SFR = complete_90+(0.4*E_stars*k_prime) - 41.1817d  ; This 41.1817 comes from Kewley et al. 2004 conv. btwn SFR and L([O II])
complete_80_SFR = complete_80+(0.4*E_stars*k_prime) - 41.1817d
complete_70_SFR = complete_70+(0.4*E_stars*k_prime) - 41.1817d
complete_60_SFR = complete_60+(0.4*E_stars*k_prime) - 41.1817d
complete_50_SFR = complete_50+(0.4*E_stars*k_prime) - 41.1817d

; I want to divide up into mass bins
lineSFR7 = dblarr(n_elements(lineSFR))
lineSFR8 = dblarr(n_elements(lineSFR))
lineSFR9 = dblarr(n_elements(lineSFR))
lineSFR10 = dblarr(n_elements(lineSFR))
z7 = dblarr(n_elements(lineSFR))
z8 = dblarr(n_elements(lineSFR))
z9 = dblarr(n_elements(lineSFR))
z10 = dblarr(n_elements(lineSFR))
for i = 0, n_elements(mass)-1 do begin
   if mass[i] ge 7.5 and mass[i] le 8.5 then begin
      lineSFR7[i] = lineSFR[i]
      z7[i] = z[i]
   endif
   if mass[i] ge 8.5 and mass[i] le 9.5 then begin
      lineSFR8[i] = lineSFR[i]
      z8[i] = z[i]
   endif
   if mass[i] ge 9.5 and mass[i] le 10.5 then begin
      lineSFR9[i] = lineSFR[i]
      z9[i] = z[i]
   endif
   if mass[i] ge 10.5 and mass[i] le 11.8 then begin
      lineSFR10[i] = lineSFR[i]
      z10[i] = z[i]
   endif
endfor

x = where(lineSFR7 eq 0)
remove, x, z7, lineSFR7
x = where(lineSFR8 eq 0)
remove, x, z8, lineSFR8
x = where(lineSFR9 eq 0)
remove, x, z9, lineSFR9
x = where(lineSFR10 eq 0)
remove, x, z10, lineSFR10

openpps, 'SFR_plots/HPS_o2sed_IRAC/lineSFRvsz' 
plot, complete_z, complete_80_SFR, xtitle = textoidl('Redshift'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), ystyle = 1, yrange = [-2, 2.75], xstyle = 1, xrange = [0.02, 0.58]
oplot, z7, lineSFR7, psym = sym(1), symsize = 1.2, color = 242
oplot, z8, lineSFR8, psym = sym(1), symsize = 1.2, color = 100
oplot, z9, lineSFR9, psym = sym(1), symsize = 1.2, color = 140
oplot, z10, lineSFR10, psym = sym(1), symsize = 1.2, color = 55
oplot, xray_z[[1, 3, 5]], xray_lineSFR[[1, 3, 5]], psym = 7, symsize = 2, color = 242
oplot, xray_z[[2, 4, 8]], xray_lineSFR[[2, 4, 8]], psym = 7, symsize = 2, color = 100
oplot, [xray_z[0]], [xray_lineSFR[0]], psym = 7, symsize = 2, color = 140
oplot, xray_z[[6, 7, 9]], xray_lineSFR[[6, 7, 9]], psym = 7, symsize = 2, color = 55
legend, ['log M < 8.5', '8.5 < log M < 9.5', '9.5 < log M < 10.5', 'log M > 10.5'], psym = [sym(1), sym(1), sym(1), sym(1)], color = [242, 100, 140, 55]
oplot, [0.53], [-1.5], psym = 3
oploterror, [0.53], [-2.25], [0], [alog10(median(10^lineSFR_up))], /hibar
oploterror, [0.53], [-2.25], [0], [alog10(median(10^lineSFR_low))], /lobar
closepps




bit1 = (EW)[where(mass lt 8.25)]
bit2 = (EW)[where(mass ge 8.25 and mass lt 8.8)]
bit3 = (EW)[where(mass ge 8.8 and mass lt 9.35)]
bit4 = (EW)[where(mass ge 9.35 and mass lt 9.8)]
bit5 = (EW)[where(mass ge 9.8 and mass lt 10.25)]
bit6 = (EW)[where(mass ge 10.25 and mass lt 10.8)]
bit7 = (EW)[where(mass ge 10.8 and mass lt 11.8)]
;bit8 = (EW)[where(mass ge 10.5 and mass lt 11.8)]
;bit9 = (EW)[where(mass ge 11.3 and mass lt 11.9)]

bit1 = percentiles(bit1, value = [0.16, 0.5, 0.84])
bit2 = percentiles(bit2, value = [0.16, 0.5, 0.84])
bit3 = percentiles(bit3, value = [0.16, 0.5, 0.84])
bit4 = percentiles(bit4, value = [0.16, 0.5, 0.84])
bit5 = percentiles(bit5, value = [0.16, 0.5, 0.84])
bit6 = percentiles(bit6, value = [0.16, 0.5, 0.84])
bit7 = percentiles(bit7, value = [0.16, 0.5, 0.84])
;bit8 = percentiles(bit8, value = [0.16, 0.5, 0.84])
;bit9 = percentiles(bit9, value = [0.16, 0.5, 0.84])

low = [bit1[0], bit2[0], bit3[0], bit4[0], bit5[0], bit6[0], bit7[0]];, bit8[0]];, bit9[0]]
med = [bit1[1], bit2[1], bit3[1], bit4[1], bit5[1], bit6[1], bit7[1]];, bit8[1]];, bit9[1]] 
high = [bit1[2], bit2[2], bit3[2], bit4[2], bit5[2], bit6[2], bit7[2]];, bit8[2]];, bit9[2]]
mass_ind_EW = [7.7, 8.45, 8.95, 9.45, 9.8, 10.35, 10.85]

;openpps, 'SFR_plots/HPS_o2sed_IRAC/EWvmass'
;plot, [12.2], [72], psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('[O II] EW ('+string(197B)+')'), xstyle = 1, xrange = [7, 13], ystyle = 1, yrange = [.1, 200], /ylog
;oplot, mass, EW, psym = 6, color = 100
;oplot, mass_ind, low, lines = 3
;oplot, mass_ind, med, lines = 3
;oplot, mass_ind, high, lines = 3
;oploterror, [12.2], [72], [avg(xerr_low)], [0], /lobar      
;oploterror, [12.2], [72], [avg(xerr_up)], [0], /hibar 
;oploterror, [12.2], [72], [0], [avg(EW_low)], /hibar
;oploterror, [12.2], [72], [0], [avg(EW_up)], /lobar
;closepps

openpps, 'SFR_plots/HPS_o2sed_IRAC/EWvmass_z'
plot, [12.2], [80], psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('[O II] Rest EW ('+string(197B)+')'), xstyle = 1, xrange = [7.2, 12], ystyle = 1, yrange = [.1, 200], /ylog
oplot, mass_lowz, EW_lowz, psym = sym(1), symsize = 1.2, color = 242
oplot, mass_mid1z, EW_mid1z, psym = sym(1), symsize = 1.2, color = 100
oplot, mass_mid2z, EW_mid2z, psym = sym(1), symsize = 1.2, color = 140
oplot, mass_hiz, EW_hiz, psym = sym(1), symsize = 1.2, color = 55
oplot, xray_mass[[1, 3, 5]], xray_EW[[1, 3, 5]], psym = 7, symsize = 2, color = 242
oplot, xray_mass[[2, 4, 8]], xray_EW[[2, 4, 8]], psym = 7, symsize = 2, color = 100
oplot, [xray_mass[0]], [xray_EW[0]], psym = 7, symsize = 2, color = 140
oplot, xray_mass[[6, 7, 9]], xray_EW[[6, 7, 9]], psym = 7, symsize = 2, color = 55
oplot, mass_ind_EW, low, lines = 3
oplot, mass_ind_EW, med, lines = 3
oplot, mass_ind_EW, high, lines = 3
oploterror, [11.55], [100], [median(xerr_low)], [0], /lobar      
oploterror, [11.55], [100], [median(xerr_up)], [0], /hibar 
oploterror, [11.55], [100], [0], [100*median(EW_low)/median(EW)], /hibar
oploterror, [11.55], [100], [0], [100*median(EW_up)/median(EW)], /lobar
legend, ['0.0 < z < 0.2', '0.2 < z < 0.35', '0.35 < z < 0.45', '0.45 < z < 0.57'], psym = [sym(1), sym(1), sym(1), sym(1)], color = [242, 100, 140, 45], pos = [7.35, .45]
closepps

openpps, 'SFR_plots/HPS_o2sed_IRAC/EW_hist'
multiplot, [3, 1],  myTitle = 'Number of [OII] Emitters', mxTitle = textoidl('[O II] Rest EW ('+string(197B)+')'), /square, mxtitoffset = 1.1
cgHistoplot, [EW_lowz, EW_mid1z], bin = 4, /line_fill, color = 'blu4', thick = 6, polycolor = 'blu4', orientation = [45, -45], xrange = [0, 85],  yrange = [0, 30],  ytitle = ' ', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black' , xtickname=['0',  '20', '40', '60', '80']
cgText, 50, 27.2, '0.0 < z 0.35', charsize = .9
multiplot
cgHistoplot, EW_mid2z,  bin = 4, /line_fill, color =  'blu4', thick = 6, polycolor = 'blu4', orientation = [45, -45], xrange = [0, 85],  yrange = [0, 30], YTickformat='(A1)',  ytitle = ' ', xcharsize = '.75', ycharsize = '.75', /oprobability, probcolor = 'black', xtickname=[ ' ', '20', '40',  '60', '80']
cgText, 50, 27.2, '0.35 < z 0.45', charsize = .9
multiplot
cgHistoplot, EW_hiz,  bin = 4, /line_fill, color = 'blu4', thick = 6, polycolor =  'blu4', orientation = [45, -45], xrange =[0, 85],  yrange = [0, 30],YTickformat='(A1)',  ytitle = ' ',xcharsize = '.75', ycharsize = '.75',  /oprobability, probcolor = 'black', xtickname=[ ' ', '20', '40', '60', '80']
cgText, 5, 27.2, '0.45 < z 0.57', charsize = .9
axis, yaxis = 1, ytitle = ' ', ycharsize = '.75', ytickname =  ['0','0.2', '0.4', '0.6', '0.8', '1.0'], yminor = 1, yticks = 5
closepps

; Want a plot where I have median mass bins and stdev for bins within
; each redshift bin
;bit1 = (lineSFR_lowz)[where(mass_lowz lt 8.5)]
;bit2 = (lineSFR_lowz)[where(mass_lowz ge 8.5 and mass_lowz lt 9)]
;bit3 = (lineSFR_lowz)[where(mass_lowz ge 9 and mass_lowz lt 10)]
;bit4 = (lineSFR_lowz)[where(mass_lowz ge 10 and mass_lowz lt 13)]
;std1 = stddev(bit1)
;std2 = stddev(bit2)
;std3 = stddev(bit3)
;std4 = stddev(bit4)
;med1 = median(bit1, /even)
;med2 = median(bit2, /even)
;med3 = median(bit3, /even)
;med4 = median(bit4, /even)
;med_lowz = [med1, med2, med3, med4]
;std_up_lowz = [med1+std1, med2+std2, med3+std3, med4+std4]
;std_low_lowz = [med1-std1, med2-std2, med3-std3, med4-std4]

;bit1 = (lineSFR_mid1z)[where(mass_mid1z lt 8.5)]
;bit2 = (lineSFR_mid1z)[where(mass_mid1z ge 8.5 and mass_mid1z lt 9)]
;bit3 = (lineSFR_mid1z)[where(mass_mid1z ge 9 and mass_mid1z lt 10)]
;bit4 = (lineSFR_mid1z)[where(mass_mid1z ge 10 and mass_mid1z lt 13)]
;std1 = stddev(bit1)
;std2 = stddev(bit2)
;std3 = stddev(bit3)
;std4 = stddev(bit4)
;med1 = median(bit1, /even)
;med2 = median(bit2, /even)
;med3 = median(bit3, /even)
;med4 = median(bit4, /even)
;med_mid1z = [med1, med2, med3, med4]
;std_up_mid1z = [med1+std1, med2+std2, med3+std3, med4+std4]
;std_low_mid1z = [med1-std1, med2-std2, med3-std3, med4-std4]

;bit1 = (lineSFR_mid2z)[where(mass_mid2z lt 8.5)]
;bit2 = (lineSFR_mid2z)[where(mass_mid2z ge 8.5 and mass_mid2z lt 9)]
;bit3 = (lineSFR_mid2z)[where(mass_mid2z ge 9 and mass_mid2z lt 10)]
;bit4 = (lineSFR_mid2z)[where(mass_mid2z ge 10 and mass_mid2z lt 13)]
;std1 = stddev(bit1)
;std2 = stddev(bit2)
;std3 = stddev(bit3)
;std4 = stddev(bit4)
;med1 = median(bit1, /even)
;med2 = median(bit2, /even)
;med3 = median(bit3, /even)
;med4 = median(bit4, /even)
;med_mid2z = [med1, med2, med3, med4]
;std_up_mid2z = [med1+std1, med2+std2, med3+std3, med4+std4]
;std_low_mid2z = [med1-std1, med2-std2, med3-std3, med4-std4]

;bit1 = (lineSFR_hiz)[where(mass_hiz lt 8.5)]
;bit2 = (lineSFR_hiz)[where(mass_hiz ge 8.5 and mass_hiz lt 9)]
;bit3 = (lineSFR_hiz)[where(mass_hiz ge 9 and mass_hiz lt 10)]
;bit4 = (lineSFR_hiz)[where(mass_hiz ge 10 and mass_hiz lt 13)]
;std1 = stddev(bit1)
;std2 = stddev(bit2)
;std3 = stddev(bit3)
;std4 = stddev(bit4)
;med1 = median(bit1, /even)
;med2 = median(bit2, /even)
;med3 = median(bit3, /even)
;med4 = median(bit4, /even)
;med_hiz = [med1, med2, med3, med4]
;std_up_hiz = [med1+std1, med2+std2, med3+std3, med4+std4]
;std_low_hiz = [med1-std1, med2-std2, med3-std3, med4-std4]
;mass_ind = [7.5, 8.5, 9.5, 11]

;loadct, 17
;openpps, 'SFR_plots/HPS_o2sed_IRAC/SFRvsM_median'
;plot, mass_ind, med_lowz,  psym = 3, xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), xstyle = 1, xrange = [7, 12.5], ystyle = 1, yrange = [-2, 3.5]
;oplot, mass_ind, med_lowz, psym = -6, color = 242
;oplot, mass_ind, med_mid1z+0.5, psym = -6, color = 100
;oplot, mass_ind, med_mid2z+1, psym = -6, color = 140
;oplot, mass_ind, med_hiz+1.5, psym = -6, color = 55
;closepps

;m = findgen(1300)/100
;y1 = 0.18*m - 1.18      ; my line fit
;y2 = 0.67*m - 6.19      ; Noeske line fit
;y3 = 0.5547*m - 5.3126  ; Lara-Lopez line fit
;openpps, 'SFR_plots/HPS_o2sed_IRAC/slopes'
;plot, m, y2, xstyle = 1, xrange = [7,12], xtitle = textoidl('log(M_* [M_{sun}])'), ytitle = textoidl('log(SFR [M_{sun} yr^{-1}])'), lines=4
;oplot, m, y3, lines = 2
;oplot, m, y1
;legend, ['This survey', 'AEGIS', 'SDSS/GAMA'], lines=[0, 4, 2]
;closepps

; Some numbers of interst
lineSSFR = lineSFR-mass
mass_std = stdev(mass)
SFR_std = stdev(lineSFR)
SSFR_std = stdev(lineSSFR)
print, 'Median mass:', median(mass[sort(mass)], /even)
print, 'Standard deviation of mass:', mass_std
print, 'Median SFR:', median(lineSFR[sort(lineSFR)], /even)
print, 'Standard deviation of SFR:', SFR_std
print, 'Median sSFR:', median(lineSSFR[sort(lineSSFR)], /even)
print, 'Standard deviation of sSFR:', SSFR_std
print, 'Mass percentiles:', percentiles(mass, value = [0.16, 0.5, 0.84])
print, 'SFR percentiles:', percentiles(lineSFR, value = [0.16, 0.5, 0.84])
print, 'sSFR percentiles:', percentiles(lineSSFR, value = [0.16, 0.5, 0.84])
print, ' '
print, 'Median mass low, mid1z:', median(([mass_lowz, mass_mid1z])[sort([mass_lowz, mass_mid1z])], /even)
print, 'Median mass mid2z:', median(mass_mid2z[sort(mass_mid2z)], /even)
print, 'Median mass hiz:', median(mass_hiz[sort(mass_hiz)], /even)
print, 'Median SFR low, mid1z:', median(([lineSFR_lowz, lineSFR_mid1z])[sort([lineSFR_lowz, lineSFR_mid1z])], /even)
print, 'Median SFR mid2z:', median(lineSFR_mid2z[sort(lineSFR_mid2z)], /even)
print, 'Median SFR hiz:', median(lineSFR_hiz[sort(lineSFR_hiz)], /even)
print, 'Median sSFR low, mid1z:', median(([lineSFR_lowz-mass_lowz, lineSFR_mid1z-mass_mid1z])[sort([lineSFR_lowz-mass_lowz, lineSFR_mid1z-mass_mid1z])], /even)
print, 'Median sSFR mid2z:', median((lineSFR_mid2z-mass_mid2z)[sort(lineSFR_mid2z-mass_mid2z)], /even)
print, 'Median sSFR hiz:', median((lineSFR_hiz-mass_hiz)[sort(lineSFR_hiz-mass_hiz)], /even)
print, ' '
print, 'Mass low, mid1z percentiles:', percentiles([mass_lowz, mass_mid1z], value = [0.16, 0.5, 0.84])
print, 'Mass mid2z percentiles:', percentiles(mass_mid2z, value = [0.16, 0.5, 0.84])
print, 'Mass hiz percentiles:', percentiles(mass_hiz, value = [0.16, 0.5, 0.84])
print, 'SFR low, mid1z percentiles:', percentiles([lineSFR_lowz, lineSFR_mid1z], value = [0.16, 0.5, 0.84])
print, 'SFR mid2z percentiles:', percentiles(lineSFR_mid2z, value = [0.16, 0.5, 0.84])
print, 'SFR hiz percentiles:', percentiles(lineSFR_hiz, value = [0.16, 0.5, 0.84])
print, 'sSFR low, mid1z percentiles:', percentiles([lineSFR_lowz-mass_lowz, lineSFR_mid1z-mass_mid1z], value = [0.16, 0.5, 0.84])
print, 'sSFR mid2z percentiles:', percentiles(lineSFR_mid2z-mass_mid2z, value = [0.16, 0.5, 0.84])
print, 'sSFR hiz percentiles:', percentiles(lineSFR_hiz-mass_hiz, value = [0.16, 0.5, 0.84])



openw, lun, 'morphology/mass_SFR.txt', /get_lun
for i = 0, n_elements(mass)-1 do begin
   printf, lun, string(IDs[i])+'   '+string(mass[i], format='(F0.3)')+'   '+string(lineSFR[i], format='(F0.4)')+'   '+string(E_BV[i], format='(F0.3)')+'   '+string(z[i], format='(F0.3)')
endfor
Free_lun, lun

;openw, lun, 'OII_stuff.txt'
;for i = 0, n_elements(age)-1 do begin
;   printf, lun, string(IDs[i])+' '+string(age[i])+' '+string(E_BV[i])
;endfor
;Free_lun, lun

;openw, lun, 'OII_coords_munics.txt'
;for i = 0, n_elements(all_OII.munics.RA)-1 do begin
;   printf, lun, string(all_OII.munics.RA[i])+' '+string(all_OII.munics.dec[i])
;endfor
;Free_lun, lun

lineSFR_lowred = dblarr(n_elements(mass))
lineSFR_mid1red = dblarr(n_elements(mass))
lineSFR_mid2red = dblarr(n_elements(mass))
lineSFR_hired = dblarr(n_elements(mass))
lineSFR_robin_lowred = dblarr(n_elements(mass))
lineSFR_robin_mid1red = dblarr(n_elements(mass))
lineSFR_robin_mid2red = dblarr(n_elements(mass))
lineSFR_robin_hired = dblarr(n_elements(mass))
UVSFR_robin_lowred = dblarr(n_elements(mass))
UVSFR_robin_mid1red = dblarr(n_elements(mass))
UVSFR_robin_mid2red = dblarr(n_elements(mass))
UVSFR_robin_hired = dblarr(n_elements(mass))

for i = 0, n_elements(z)-1 do begin
   if E_BV[i] ge 0 and E_BV[i] lt 0.25 then begin
      lineSFR_lowred[i] = lineSFR[i]
      lineSFR_robin_lowred[i] = lineSFR_robin[i]
      UVSFR_robin_lowred[i] = UVSFR_robin[i]
   endif
   if E_BV[i] ge 0.25 and E_BV[i] lt 0.5 then begin
      lineSFR_mid1red[i] = lineSFR[i]
      lineSFR_robin_mid1red[i] = lineSFR_robin[i]
      UVSFR_robin_mid1red[i] = UVSFR_robin[i]
   endif
   if E_BV[i] ge 0.5 and E_BV[i] lt 0.75 then begin
      lineSFR_mid2red[i] = lineSFR[i]
      lineSFR_robin_mid2red[i] = lineSFR_robin[i]
      UVSFR_robin_mid2red[i] = UVSFR_robin[i]
   endif
   if E_BV[i] ge 0.75 and E_BV[i] le 1.0 then begin
      lineSFR_hired[i] = lineSFR[i]
      lineSFR_robin_hired[i] = lineSFR_robin[i]
      UVSFR_robin_hired[i] = UVSFR_robin[i]
   endif
endfor





stop
end
