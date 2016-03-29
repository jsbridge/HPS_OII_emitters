;+
;
; Tabulates chi squared values for the SED fits
;
;-
function chi_squared

restore, 'ascii_templates/chi_template.sav'

;cos_file = findfile('SED_fitting/cos_80k_like/*.likestats')
;cos_file = findfile('SED_fitting/cosmos_IRAC_like/*.likestats')
;goods_file = findfile('SED_fitting/goods_80k_like/*.likestats')
;goods_file = findfile('SED_fitting/goods_final_like/*.likestats')
;munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
;xmm_file = findfile('SED_fitting/xmm_80k_like/*.likestats')
;xmm_file = findfile('SED_fitting/xmm_IRAC_like/*.likestats')

cos_file = findfile('SED_fitting/cosmos_IRAC_like/*.likestats')
goods_file = findfile('SED_fitting/goods_final_like/*.likestats')
munics_file = findfile('SED_fitting/munics_80k_like/*.likestats')
xmm_file = findfile('SED_fitting/xmm_IRAC_like/*.likestats')

c_chi = dblarr(n_elements(cos_file))
g_chi = dblarr(n_elements(goods_file))
m_chi = dblarr(n_elements(munics_file))
x_chi = dblarr(n_elements(xmm_file))

for i = 0, n_elements(cos_file)-1 do begin
   c_data = read_ascii(cos_file[i], template = chi_template) 
   c_chi[i] = c_data.field6[0]  
endfor
for i = 0, n_elements(goods_file)-1 do begin
   g_data = read_ascii(goods_file[i], template = chi_template) 
   g_chi[i] = g_data.field6[0]  
endfor
for i = 0, n_elements(munics_file)-1 do begin
   m_data = read_ascii(munics_file[i], template = chi_template) 
   m_chi[i] = m_data.field6[0]  
endfor
for i = 0, n_elements(xmm_file)-1 do begin
   x_data = read_ascii(xmm_file[i], template = chi_template) 
   x_chi[i] = x_data.field6[0] 
endfor

chi = 2d * [x_chi, m_chi, c_chi, g_chi]

return, chi

end
