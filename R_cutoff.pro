;+
;
; This code implements an R-1 cutoff of 0.2 for the SED fits. It
; returns a structure with arrays for each field that contain the IDs
; of the objects that fit the cutoff.  Each array is the length of the
; full field of OII emitters (ie, cosmos array is 110 elements), but
; the objects that don't fit the criteria have their ID numbers set to 0
;
;-
function R_cutoff

all_OII = all_OII()


; Find the .converge files for each OII emitter in HPS
restore, 'ascii_templates/converge_template.sav'
;cos_file = findfile('SED_fitting/cos_80k_conv/*.converge')
cos_file = findfile('SED_fitting/cosmos_IRAC_conv/*.converge')
;goods_file = findfile('SED_fitting/goods_80k_conv/*.converge')
goods_file = findfile('SED_fitting/goods_final_conv/*.converge')
munics_file = findfile('SED_fitting/munics_80k_conv/*.converge')
;xmm_file = findfile('SED_fitting/xmm_80k_conv/*.converge')
xmm_file = findfile('SED_fitting/xmm_IRAC_conv/*.converge')

; Read in each file and assign the R-1 value of each parameter
; (E(B-V), mass, etc) for each file.  What I want to do is not assign
; them to an array, but just read the value, check if it's
; greater than 0.2 and if not, keep track of the ID number for which
; ones satisfy the criterion and which don't
c_ID = dblarr(n_elements(cos_file))
g_ID = dblarr(n_elements(goods_file))
m_ID = dblarr(n_elements(munics_file))
x_ID = dblarr(n_elements(xmm_file))
c_IDb = dblarr(n_elements(cos_file))
g_IDb = dblarr(n_elements(goods_file))
m_IDb = dblarr(n_elements(munics_file))
x_IDb = dblarr(n_elements(xmm_file))
c_R = dblarr(n_elements(cos_file))
g_R = dblarr(n_elements(goods_file))
m_R = dblarr(n_elements(munics_file))
x_R = dblarr(n_elements(xmm_file))

for i = 0, n_elements(cos_file)-1 do begin
   c_data = read_ascii(cos_file[i], template = converge_template) 
   c_R[i] = c_data.field2[3] 
   if c_R[i] le 0.2 then begin
      c_ID[i] = 1
   endif else begin
      c_IDb[i] = 1
   endelse
endfor
for i = 0, n_elements(goods_file)-1 do begin
   g_data = read_ascii(goods_file[i], template = converge_template) 
   g_R[i] = g_data.field2[3] 
   if g_R[i] le 0.2 then begin
      g_ID[i] = 1
   endif else begin
      g_IDb[i] = 1
   endelse
endfor
for i = 0, n_elements(munics_file)-1 do begin
   m_data = read_ascii(munics_file[i], template = converge_template) 
   m_R[i] = m_data.field2[3] 
   if m_R[i] le 0.2 then begin
      m_ID[i] = 1
   endif else begin
      m_IDb[i] = 1
   endelse
endfor
for i = 0, n_elements(xmm_file)-1 do begin
   x_data = read_ascii(xmm_file[i], template = converge_template) 
   x_R[i] = x_data.field2[3] 
   if x_R[i] le 0.2 then begin
      x_ID[i] = 1
   endif else begin
      x_IDb[i] = 1
   endelse
endfor

;This little section pulls out the IDs of the ones that don't
;pass muster and prints them to the screen, just in case they need
;further investigating
a = [c_IDb*all_OII.cosmos.id, g_IDb*all_OII.goods.id, m_IDb*all_OII.munics.id, x_IDb*all_OII.xmm.id]
b = where(a eq 0)
if n_elements(a) ne n_elements(b) then begin
   remove, b, a
   print, fix(a)
endif

c_ID = c_ID * all_OII.cosmos.id
g_ID = g_ID * all_OII.goods.id
m_ID = m_ID * all_OII.munics.id
x_ID = x_ID * all_OII.xmm.id

return, create_struct('cosmos', c_ID, 'goods', g_ID, 'munics', m_ID, 'xmm', x_ID, 'c_R', c_R, 'g_R', g_R, 'm_R', m_R, 'x_R', x_R)

end
