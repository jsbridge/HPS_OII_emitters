;+
;
; The purpose of this function is to recalculate the SFRs from the
; line luminosities using the E(B-V) from the SED fitting.  These will
; be used in place of Robin's SFRs
;
; From Calzetti: F_i = F_o * 10^(0.4 * E_s * k_prime)
;
;-
function luminosity

all_OII = all_OII()

restore, 'ascii_templates/likestats_template.sav'

;cos_file = findfile('SED_fitting/cos_80k_const/*.likestats')
;goods_file = findfile('SED_fitting/goods_80k_const/*.likestats')
;munics_file = findfile('SED_fitting/munics_80k_const/*.likestats')
;xmm_file = findfile('SED_fitting/xmm_80k_const/*.likestats')

;cos_file = findfile('SED_fitting/cos_80k_like/*.likestats')
;goods_file = findfile('SED_fitting/goods_80k_like/*.likestats')
;munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
;xmm_file = findfile('SED_fitting/xmm_80k_like/*.likestats')

cos_file = findfile('SED_fitting/cosmos_IRAC_like/*.likestats')
goods_file = findfile('SED_fitting/goods_final_like/*.likestats')
munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
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
c_E_BV = dblarr(n_elements(cos_file))
g_E_BV = dblarr(n_elements(goods_file))
m_E_BV = dblarr(n_elements(munics_file))
x_E_BV = dblarr(n_elements(xmm_file))
c_E_BV_low = dblarr(n_elements(cos_file))
c_E_BV_up = dblarr(n_elements(cos_file))
g_E_BV_low = dblarr(n_elements(goods_file))
g_E_BV_up = dblarr(n_elements(goods_file))
m_E_BV_low = dblarr(n_elements(munics_file))
m_E_BV_up = dblarr(n_elements(munics_file))
x_E_BV_low = dblarr(n_elements(xmm_file))
x_E_BV_up = dblarr(n_elements(xmm_file))

for i = 0, n_elements(cos_file)-1 do begin
   c_data = read_ascii(cos_file[i], template = likestats_template) 
   c_E_BV[i] = c_data.field2[2]  
   c_E_BV_low[i] = c_data.field3[2]
   c_E_BV_up[i] = c_data.field4[2]
   c_mass[i] = c_data.field2[3] 
   c_mass_up[i] = c_data.field4[3]-c_data.field2[3]  ; 1 sigma
   c_mass_low[i] = c_data.field2[3]-c_data.field3[3]  ; 1 sigma
endfor
for i = 0, n_elements(goods_file)-1 do begin
   g_data = read_ascii(goods_file[i], template = likestats_template) 
   g_E_BV[i] = g_data.field2[2]  
   g_E_BV_low[i] = g_data.field3[2]
   g_E_BV_up[i] = g_data.field4[2]
   g_mass[i] = g_data.field2[3] 
   g_mass_up[i] = g_data.field4[3]-g_data.field2[3]  ; 1 sigma
   g_mass_low[i] = g_data.field2[3]-g_data.field3[3]  ; 1 sigma
endfor
for i = 0, n_elements(munics_file)-1 do begin
   m_data = read_ascii(munics_file[i], template = likestats_template) 
   m_E_BV[i] = m_data.field2[2]  
   m_E_BV_low[i] = m_data.field3[2]
   m_E_BV_up[i] = m_data.field4[2]
   m_mass[i] = m_data.field2[3]  
   m_mass_up[i] = m_data.field4[3]-m_data.field2[3]  ; 1 sigma
   m_mass_low[i] = m_data.field2[3]-m_data.field3[3]  ; 1 sigma
endfor
for i = 0, n_elements(xmm_file)-1 do begin
   x_data = read_ascii(xmm_file[i], template = likestats_template) 
   x_E_BV[i] = x_data.field2[2] 
   x_E_BV_low[i] = x_data.field3[2]
   x_E_BV_up[i] = x_data.field4[2] 
   x_mass[i] = x_data.field2[3]  
   x_mass_up[i] = x_data.field4[3]-x_data.field2[3]  ; 1 sigma
   x_mass_low[i] = x_data.field2[3]-x_data.field3[3]  ; 1 sigma
endfor

E_BV_sed =  [x_E_BV, m_E_BV, c_E_BV, g_E_BV]
E_BV_upsed = [x_E_BV_up, m_E_BV_up, c_E_BV_up, g_E_BV_up] - E_BV_sed
E_BV_lowsed = E_BV_sed - [x_E_BV_low, m_E_BV_low, c_E_BV_low, g_E_BV_low]
mass = [x_mass, m_mass, c_mass, g_mass]
mass = mass*alog10(2.71828182846)
mass_up = [x_mass_up, m_mass_up, c_mass_up, g_mass_up]*alog10(2.71828182846)
mass_low = [x_mass_low, m_mass_low, c_mass_low, g_mass_low]*alog10(2.71828182846)

; Read in HPStable3 fluxes
restore, 'ascii_templates/table3_template.sav'
HPS = read_ascii('data/HPStable3.txt', template = table3_template)

; Picks out OII emitters
readcol, 'data/robin_sfr.dat', robin_OIIid, robin_EBV, robin_OII_SFR, format = 'I,X,X,D,X,X,D'
ind = intarr(n_elements(robin_OIIid))
for i = 0, n_elements(robin_OIIid)-1 do begin
   ind[i] = where(robin_OIIid[i] eq HPS.field01)
endfor

; So, I have a bunch of otions with which to fill in the missing
; E(B-V) values for which I don't have defined UV slopes
; One is to use the reddening from Josh's SED fitting, or I
; could use my SED fitting, which I don't think is a good idea
; (I tried it), or I could use the Garns and Best relation,
; which I'm going to go with, I think.

;E_BV_sed = robin_EBV
E_BV_sed = garn_best(mass)

; This whole section basically forms E_BV with the 215 UV slope
; values, and then populates the rest of the array with the SED/G&B values
readcol, 'OII_uv_slopes.dat', obj_new, E_BV_new, E_BV_newerr, format = 'I,X,X,D,D'

x = intarr(n_elements(E_BV_new))
for i = 0, n_elements(E_BV_new)-1 do begin
   x[i] = where(obj_new[i] eq robin_OIIid, complement = y)
endfor

E_BV = fltarr(n_elements(E_BV_sed))
E_BV_up = fltarr(n_elements(E_BV_sed))
E_BV_low = fltarr(n_elements(E_BV_sed))
for i = 0, n_elements(E_BV_new)-1 do begin
   E_BV[x[i]] = E_BV_new[i]
   E_BV_up[x[i]] = E_BV_newerr[i]
   E_BV_low[x[i]] = E_BV_newerr[i]
endfor

p = where(E_BV eq 0)
s = where(mass lt 8.5)

r = intarr(n_elements(s))
for i = 0, n_elements(s)-1 do begin
   r[i] = where(s[i] eq p)
endfor

b = where(r eq -1)
remove, b, r
c = p[r]

E_BV[c] = min(E_BV_sed)

for i = 0, n_elements(E_BV_sed)-1 do begin
  if E_BV[i] eq 0 then begin
         E_BV[i] = E_BV_sed[i]
         E_BV_up[i] = E_BV_upsed[i]
         E_BV_low[i] = E_BV_lowsed[i]
  endif
endfor


OII_flux = HPS.field11[ind]*10^(-17d)  ; mW/m^2  - No conversion needed to get this in ergs/s/cm^2
OII_flux_up = HPS.field13[ind]*10^(-17d)
OII_flux_low = HPS.field12[ind]*10^(-17d)

;for i = 0, n_elements(OII_flux)-1 do begin
;   print, OII_flux[i]*10^18d
;   print, OII_flux_up[i]*10^18d
;   print, OII_flux_low[i]*10^18d
;endfor

;From Robin's new code calzetti.f
wp = 0.3727
R_OII = 1.78 + (1.17 * ((0.011 / (wp^3)) - (0.198 / (wp^2)) + (1.509 / wp) - 2.156)) 
A_OII = R_OII/0.44    ; As per Calzetti for emission line extinction

;From Calzetti et al.
k_prime = 2.659d*((1.509d/0.3727) - (0.198d/0.3727^2) + (0.011d/0.3727^3) - 2.156d) + 4.05d
k_prime_err = 0.8d

k_prime = A_OII

; Progpagate the error - this is the variance of E(B-V)*k'
E_k_prime = E_BV*k_prime
E_k_prime_up = E_BV*k_prime * sqrt((E_BV_up/E_BV)^2 + (k_prime_err/k_prime)^2)
E_k_prime_low = E_BV*k_prime * sqrt((E_BV_low/E_BV)^2 + (k_prime_err/k_prime)^2)

; Calculate the 10^Esk' part, and propagate the error (using
; approx. from wikipedia)
log_part = 10^(0.4d*E_k_prime)
log_part_low = log_part * alog(10) * 0.4d*E_k_prime_low
log_part_up = log_part * alog(10) * 0.4d*E_k_prime_up

; Now, include the F_o (flux, in ergs/s/cm^2) from Calzetti,
; propagating the error from the measurement and the error above

dered_flux = OII_flux * log_part
dered_flux_low = dered_flux * sqrt((OII_flux_low/OII_flux)^2 + (log_part_low/log_part)^2)
dered_flux_up = dered_flux * sqrt((OII_flux_up/OII_flux)^2 + (log_part_up/log_part)^2)

; Now need to convert from flux to luminosity where R comes from the
; redshift...R is luminosity distance, so need lumdist.pro!

all_OII = all_OII()
z = [all_OII.xmm.z, all_OII.munics.z, all_OII.cosmos.z, all_OII.goods.z]
OII_lum = dblarr(n_elements(OII_flux))
OII_lum_up = dblarr(n_elements(OII_flux_up))
OII_lum_low = dblarr(n_elements(OII_flux_low))
OII_lum_red = dblarr(n_elements(OII_flux))

for i = 0, n_elements(OII_flux)-1 do begin
   d = lumdist(z[i], /silent)*3.085677581d24     ; lumdist gives in Mpc, so need to convert that to cm
   OII_lum[i] = dered_flux[i]*4*!pi*d^2
   OII_lum_red[i] = OII_flux[i]*4*!pi*d^2
   OII_lum_up[i] = dered_flux_up[i]*4*!pi*d^2
   OII_lum_low[i] = dered_flux_low[i]*4*!pi*d^2
endfor

SFR_red =  6.58d-42*OII_lum_red   ; M_sun/yr
SFR = 6.58d-42*OII_lum   ; M_sun/yr
SFR_up = SFR * sqrt((1.65d-42/6.58d-42)^2 + (OII_lum_up/OII_lum)^2)
SFR_low = SFR * sqrt((1.65d-42/6.58d-42)^2 + (OII_lum_low/OII_lum)^2)



return, create_struct('robin_nums', robin_OIIid, 'nums', HPS.field01[ind],  'SFR', alog10(SFR), 'SFR_up', abs(SFR_up/(SFR*alog(10))), 'SFR_low', abs(SFR_low/(SFR*alog(10))), 'E_BV', E_BV, 'lum', OII_lum, 'SFR_red', alog10(SFR_red), 'number', x)

end
