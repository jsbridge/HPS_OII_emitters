;+
;
; Returns a structures with the OII emitter data from HPS, such as HPS ID, RA, dec, and z
;
;-
function all_OII

readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

restore, 'ascii_templates/HPS4_template.sav'
HPS_data = read_ascii('data/HPStable4.txt', template=template)

restore, 'ascii_templates/table3_template.sav'
HPS = read_ascii('data/HPStable3.txt', template = table3_template)

name = HPS_data.field11
z = HPS_data.field12
id = HPS_data.field01
Rmag = HPS_data.field03
xray = HPS_data.field14
EW = HPS_data.field05
EW_low = HPS_data.field06
EW_up = HPS_data.field07
Flux = HPS.field11
;ind = where(HPSRA eq -99)
;remove, ind, HPSRA, HPSdec, z, name, id, Rmag

O_index = where(name eq '[OII]')   ; picks out all OII emitters in HPS
xray = xray[O_index]

cos_OII = O_index[where(O_index ge 138 and O_index le 330)]   ; picks out the COSMOS OII emitters from all OII emitters
goods_OII = O_index[where(O_index ge 331 and O_index le 478)]
munics_OII = O_index[where(O_index ge 41 and O_index le 137)]
xmm_OII = O_index[where(O_index ge 0 and O_index le 40)]

cos_RA = HPSRA[142:330]   ; picks out the RA of everything in the COSMOS section of HPS
cos_dec = HPSdec[142:330]  ; same as above, but dec
g_RA = HPSRA[331:478]
g_dec = HPSdec[331:478]
m_RA = HPSRA[41:141]
m_dec = HPSdec[41:141]
x_RA = HPSRA[0:40]
x_dec = HPSdec[0:40]

cos_OIIRA = cos_RA[cos_OII - 142]  ; picks out the RA of all COSMOS OII emitters in HPS
cos_OIIdec = cos_dec[cos_OII - 142]  ; same as above, but dec
g_OIIRA = g_RA[goods_OII - 331]
g_OIIdec = g_dec[goods_OII - 331]
m_OIIRA = m_RA[munics_OII - 41]
m_OIIdec = m_dec[munics_OII - 41]
x_OIIRA = x_RA[xmm_OII]
x_OIIdec = x_dec[xmm_OII]

cos_OIIz = z[cos_OII]
g_OIIz = z[goods_OII]
m_OIIz = z[munics_OII]
x_OIIz = z[xmm_OII]

cos_OIIid = id[cos_OII]
g_OIIid = id[goods_OII]
m_OIIid = id[munics_OII]
x_OIIid = id[xmm_OII]

cos_Rmag = Rmag[cos_OII]
g_Rmag = Rmag[goods_OII]
m_Rmag = Rmag[munics_OII]
x_Rmag = Rmag[xmm_OII]

cos_EW = EW[cos_OII]
g_EW = EW[goods_OII]
m_EW = EW[munics_OII]
x_EW = EW[xmm_OII]

cos_EW_low = EW_low[cos_OII]
g_EW_low = EW_low[goods_OII]
m_EW_low = EW_low[munics_OII]
x_EW_low = EW_low[xmm_OII]

cos_EW_up = EW_up[cos_OII]
g_EW_up = EW_up[goods_OII]
m_EW_up = EW_up[munics_OII]
x_EW_up = EW_up[xmm_OII]

cos_flux = Flux[cos_OII]
g_flux = Flux[goods_OII]
m_flux = Flux[munics_OII]
x_flux = Flux[xmm_OII]

ind = where(xray ne '          NaN')

cosmos_OII = create_struct('ID', cos_OIIid, 'RA', cos_OIIRA, 'dec', cos_OIIdec, 'z', cos_OIIz, 'Rmag', cos_Rmag, 'EW', cos_EW, 'EW_low', cos_EW_low, 'EW_up', cos_EW_up, 'Flux', cos_flux)
goods_OII = create_struct('ID', g_OIIid, 'RA', g_OIIRA, 'dec', g_OIIdec, 'z', g_OIIz, 'Rmag', g_Rmag, 'EW', g_EW, 'EW_low', g_EW_low, 'EW_up', g_EW_up, 'Flux', g_flux)
munics_OII = create_struct('ID', m_OIIid, 'RA', m_OIIRA, 'dec', m_OIIdec, 'z', m_OIIz, 'Rmag', m_Rmag, 'EW', m_EW, 'EW_low', m_EW_low, 'EW_up', m_EW_up, 'Flux', m_flux)
xmm_OII = create_struct('ID', x_OIIid, 'RA', x_OIIRA, 'dec', x_OIIdec, 'z', x_OIIz, 'Rmag', x_Rmag, 'EW', x_EW, 'EW_low', x_EW_low, 'EW_up', x_EW_up, 'Flux', x_flux)

OIIid = [x_OIIid, m_OIIid, cos_OIIid, g_OIIid]
xray_id = OIIid[ind]

all_OII = create_struct('cosmos', cosmos_OII, 'goods', goods_OII, 'munics', munics_OII, 'xmm', xmm_OII, 'xray_ind', ind)


;openw, lun, '../HPS_OII_coords.txt', /get_lun
;for i = 0, n_elements(all_OII.xmm.id)-1 do begin
;   printf, lun, string(all_OII.xmm.id[i])+' '+string(all_OII.xmm.RA[i])+' '+string(all_OII.xmm.dec[i])+' '+string(all_OII.xmm.z[i])+' '+string(all_OII.xmm.Flux[i])
;endfor
;for i = 0, n_elements(all_OII.munics.id)-1 do begin
;   printf, lun, string(all_OII.munics.id[i])+' '+string(all_OII.munics.RA[i])+' '+string(all_OII.munics.dec[i])+' '+string(all_OII.munics.z[i])+' '+string(all_OII.munics.Flux[i])
;endfor
;for i = 0, n_elements(all_OII.cosmos.id)-1 do begin
;   printf, lun, string(all_OII.cosmos.id[i])+' '+string(all_OII.cosmos.RA[i])+' '+string(all_OII.cosmos.dec[i])+' '+string(all_OII.cosmos.z[i])+' '+string(all_OII.cosmos.Flux[i])
;endfor
;for i = 0, n_elements(all_OII.goods.id)-1 do begin
;   printf, lun, string(all_OII.goods.id[i])+' '+string(all_OII.goods.RA[i])+' '+string(all_OII.goods.dec[i])+' '+string(all_OII.goods.z[i])+' '+string(all_OII.goods.Flux[i])
;endfor
;Free_lun, lun

return, all_OII

end
