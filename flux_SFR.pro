function flux_SFR

all_OII = all_OII()

restore, 'ascii_templates/likestats_template.sav'

cos_file = findfile('SED_fitting/cos_80k_like/*.likestats')
goods_file = findfile('SED_fitting/goods_80k_like/*.likestats')
munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
xmm_file = findfile('SED_fitting/xmm_80k_like/*.likestats')

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
endfor
for i = 0, n_elements(goods_file)-1 do begin
   g_data = read_ascii(goods_file[i], template = likestats_template) 
   g_E_BV[i] = g_data.field2[2]  
   g_E_BV_low[i] = g_data.field3[2]
   g_E_BV_up[i] = g_data.field4[2]
endfor
for i = 0, n_elements(munics_file)-1 do begin
   m_data = read_ascii(munics_file[i], template = likestats_template) 
   m_E_BV[i] = m_data.field2[2]  
   m_E_BV_low[i] = m_data.field3[2]
   m_E_BV_up[i] = m_data.field4[2]
endfor
for i = 0, n_elements(xmm_file)-1 do begin
   x_data = read_ascii(xmm_file[i], template = likestats_template) 
   x_E_BV[i] = x_data.field2[2] 
   x_E_BV_low[i] = x_data.field3[2]
   x_E_BV_up[i] = x_data.field4[2] 
endfor

E_BV =  [x_E_BV, m_E_BV, c_E_BV, g_E_BV]
E_BV_up = [x_E_BV_up, m_E_BV_up, c_E_BV_up, g_E_BV_up] - E_BV
E_BV_low = E_BV - [x_E_BV_low, m_E_BV_low, c_E_BV_low, g_E_BV_low]

; Read in HPStable3 fluxes
restore, 'ascii_templates/table3_template.sav'
HPS = read_ascii('data/HPStable3.txt', template = table3_template)

; Picks out OII emitters
readcol, 'data/robin_sfr.dat', robin_OIIid, robin_EBV, robin_OII_SFR, format = 'I,X,X,D,X,X,D'
ind = intarr(n_elements(robin_OIIid))
for i = 0, n_elements(robin_OIIid)-1 do begin
   ind[i] = where(robin_OIIid[i] eq HPS.field01)
endfor

;E_BV = robin_EBV

log_OII_flux = alog10((HPS.field11[ind])*10^(-17d))  ; mW/m^2  - No conversion needed to get this in ergs/s/cm^2
log_OII_flux_up = alog10((HPS.field13[ind])*10^(-17d))
log_OII_flux_low = alog10((HPS.field12[ind])*10^(-17d))

;From Robin's new code calzetti.f
wp = 0.3727d
R_OII = 1.78d + (1.17d * ((0.011d / (wp^3)) - (0.198d / (wp^2)) + (1.509d / wp) - 2.156d)) 
A_OII = R_OII*E_BV/0.44
R_OII_err = 0.8d

; Progpagate the error
A_OII_up = A_OII * sqrt((E_BV_up/E_BV)^2 + (R_OII_err/R_OII)^2)
A_OII_low = A_OII * sqrt((E_BV_low/E_BV)^2 + (R_OII_err/R_OII)^2)


OII_log_flux = log_OII_flux + (A_OII/2.5d)
OII_log_flux_up = sqrt((log_OII_flux_up)^2 + (A_OII_up/2.5d)^2)
OII_log_flux_low = sqrt((log_OII_flux_low)^2 + (A_OII_low/2.5d)^2)


; Now need to convert from flux to luminosity where R comes from the
; redshift...R is luminosity distance, so need lumdist.pro!

all_OII = all_OII()
z = [all_OII.xmm.z, all_OII.munics.z, all_OII.cosmos.z, all_OII.goods.z]
log_OII_lum = dblarr(n_elements(OII_log_flux))
log_OII_lum_up = dblarr(n_elements(OII_log_flux_up))
log_OII_lum_low = dblarr(n_elements(OII_log_flux_low))
d = dblarr(n_elements(z))

for i = 0, n_elements(OII_log_flux)-1 do begin
   d[i] = lumdist(z[i], /silent)            ; lumdist gives in Mpc, so need to convert that to cm
   log_OII_lum[i] = 50.0779d + (2d * alog10(d[i])) + OII_log_flux[i]
   log_OII_lum_up[i] = OII_log_flux_up[i]
   log_OII_lum_low[i] = OII_log_flux_low[i]
endfor


SFR = log_OII_lum - 41.1818   ; M_sun/yr
SFR_up = sqrt((-41.7825)^2 + (log_OII_lum_up)^2)
SFR_low = sqrt((-41.7825)^2 + (log_OII_lum_low)^2)

return, create_struct('nums', robin_OIIid, 'SFR', SFR, 'SFR_up', SFR_up, 'SFR_low', SFR_low, 'E_BV', E_BV)


end
