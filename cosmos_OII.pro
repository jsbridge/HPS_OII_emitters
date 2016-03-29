;+
;
; Returns a structure with the COSMOS OII emitter data from HPS, such as HPS ID, RA, dec, and z
;
;-
function cosmos_OII

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
cos_OII = O_index[where(O_index ge 140 and O_index le 323)]   ; picks out the COSMOS OII emitters from all OII emitters
cos_RA = HPSRA[140:323]   ; picks out the RA of everything in the COSMOS section of HPS
cos_dec = HPSdec[140:323]  ; same as above, but dec
cos_OIIRA = cos_RA[cos_OII - 140]  ; picks out the RA of all COSMOS OII emitters in HPS
cos_OIIdec = cos_dec[cos_OII - 140]  ; same as above, but dec
cos_OIIz = z[cos_OII]
cos_OIIid = id[cos_OII]
cos_Rmag = Rmag[cos_OII]

cosmos_OII = create_struct('ID', cos_OIIid, 'RA', cos_OIIRA, 'dec', cos_OIIdec, 'z', cos_OIIz, 'Rmag', cos_Rmag)

return, cosmos_OII

end
