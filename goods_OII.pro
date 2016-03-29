;+
;
; Returns a structure with the GOODS OII emitter data from HPS, such as HPS ID, RA, dec, and z
;
;-
function goods_OII

readcol, 'data/HPS_contcoords.txt', HPSRA, HPSdec, format = 'D,D', comment=';'

restore, 'ascii_templates/HPS4_template.sav'
HPS_data = read_ascii('data/HPStable4.txt', template=template)

oxy = HPS_data.field11
z = HPS_data.field12
id = HPS_data.field01
Rmag = HPS_data.field03
ind = where(HPSRA eq -99)
remove, ind, HPSRA, HPSdec, z, oxy, id, Rmag

O_index = where(oxy eq '[OII]')   ; picks out all OII emitters in HPS
goods_OII = O_index[where(O_index ge 324 and O_index le 467)]   ; picks out the GOODS OII emitters from all OII emitters
goods_RA = HPSRA[324:467]   ; picks out the RA of everything in the GOODS section of HPS
goods_dec = HPSdec[324:467]  ; same as above, but dec
goods_OIIRA = goods_RA[goods_OII - 324]  ; picks out the RA of all GOODS OII emitters in HPS
goods_OIIdec = goods_dec[goods_OII - 324]  ; same as above, but dec
goods_OIIz = z[goods_OII]
goods_OIIid = id[goods_OII]
goods_Rmag = Rmag[goods_OII]

goods_OII = create_struct('ID', goods_OIIid, 'RA', goods_OIIRA, 'dec', goods_OIIdec, 'z', goods_OIIz, 'Rmag', goods_Rmag)

return, goods_OII

end
