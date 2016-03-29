function pirzkal

readcol, 'data/PN_Mass_SFR/global_SFR_hopkins/global_SFR_mass_0.2_0.5.txt', P_mass1, P_SFR1, P_z1, format = 'D, D, D'

Pbit1 = (P_SFR1)[where(P_mass1 lt 6.7)]
Pbit2 = (P_SFR1)[where(P_mass1 ge 6.7 and P_mass1 lt 7.1)]
Pbit3 = (P_SFR1)[where(P_mass1 ge 7.1 and P_mass1 lt 7.6)]
Pbit4 = (P_SFR1)[where(P_mass1 ge 7.6 and P_mass1 lt 8.1)]
Pbit5 = (P_SFR1)[where(P_mass1 ge 8.1 and P_mass1 lt 8.6)]
Pbit6 = (P_SFR1)[where(P_mass1 ge 8.6 and P_mass1 lt 9.1)]
Pbit7 = (P_SFR1)[where(P_mass1 ge 9.1 and P_mass1 lt 9.6)]
Pbit8 = (P_SFR1)[where(P_mass1 ge 9.6 and P_mass1 lt 10.1)]
Pbit9 = (P_SFR1)[where(P_mass1 ge 10.1 and P_mass1 lt 10.6)]
Pbit10 = (P_SFR1)[where(P_mass1 ge 10.6 and P_mass1 lt 11.1)]

Pstd1 = stddev(Pbit1)
Pstd2 = stddev(Pbit2)
Pstd3 = stddev(Pbit3)
Pstd4 = stddev(Pbit4)
Pstd5 = stddev(Pbit5)
Pstd6 = stddev(Pbit6)
Pstd7 = stddev(Pbit7)
Pstd8 = stddev(Pbit8)
Pstd9 = stddev(Pbit9)
Pstd10 = stddev(Pbit10)

Pmed1 = median(Pbit1, /even)
Pmed2 = median(Pbit2, /even)
Pmed3 = median(Pbit3, /even)
Pmed4 = median(Pbit4, /even)
Pmed5 = median(Pbit5, /even)
Pmed6 = median(Pbit6, /even)
Pmed7 = median(Pbit7, /even)
Pmed8 = median(Pbit8, /even)
Pmed9 = median(Pbit9, /even)
Pmed10 = median(Pbit10, /even)

Pmed = [Pmed1, Pmed2, Pmed3, Pmed4, Pmed5, Pmed6, Pmed7, Pmed8, Pmed9, Pmed10]
Pstd_up = [Pmed1+Pstd1, Pmed2+Pstd2, Pmed3+Pstd3, Pmed4+Pstd4, Pmed5+Pstd5, Pmed6+Pstd6, Pmed7+Pstd7, Pmed8+Pstd8, Pmed9+Pstd9, Pmed10+Pstd10]
Pstd_low = [Pmed1-Pstd1, Pmed2-Pstd2, Pmed3-Pstd3, Pmed4-Pstd4, Pmed5-Pstd5, Pmed6-Pstd6, Pmed7-Pstd7, Pmed8-Pstd8, Pmed9-Pstd9, Pmed10-Pstd10]
Pmass_ind = [6.1, 6.9, 7.35, 7.85, 8.35, 8.85, 9.35, 9.85, 10.35, 10.85]

one_Pmed = Pmed
one_Pstd_up = Pstd_up
one_Pstd_low = Pstd_low
one_Pmass_ind = Pmass_ind

readcol, 'data/PN_Mass_SFR/global_SFR_hopkins/global_SFR_mass_0.5_0.7.txt', P_mass1, P_SFR1, P_z1, format = 'D, D, D'

Pbit1 = (P_SFR1)[where(P_mass1 lt 7.2)]
Pbit2 = (P_SFR1)[where(P_mass1 ge 7.2 and P_mass1 lt 7.7)]
Pbit3 = (P_SFR1)[where(P_mass1 ge 7.7 and P_mass1 lt 8.2)]
Pbit4 = (P_SFR1)[where(P_mass1 ge 8.2 and P_mass1 lt 8.7)]
Pbit5 = (P_SFR1)[where(P_mass1 ge 8.7 and P_mass1 lt 9.2)]
Pbit6 = (P_SFR1)[where(P_mass1 ge 9.2 and P_mass1 lt 9.7)]
Pbit7 = (P_SFR1)[where(P_mass1 ge 9.7 and P_mass1 lt 10.2)]
Pbit8 = (P_SFR1)[where(P_mass1 ge 10.2 and P_mass1 lt 11.7)]

Pstd1 = stddev(Pbit1)
Pstd2 = stddev(Pbit2)
Pstd3 = stddev(Pbit3)
Pstd4 = stddev(Pbit4)
Pstd5 = stddev(Pbit5)
Pstd6 = stddev(Pbit6)
Pstd7 = stddev(Pbit7)
Pstd8 = stddev(Pbit8)

Pmed1 = median(Pbit1, /even)
Pmed2 = median(Pbit2, /even)
Pmed3 = median(Pbit3, /even)
Pmed4 = median(Pbit4, /even)
Pmed5 = median(Pbit5, /even)
Pmed6 = median(Pbit6, /even)
Pmed7 = median(Pbit7, /even)
Pmed8 = median(Pbit8, /even)

Pmed = [Pmed1, Pmed2, Pmed3, Pmed4, Pmed5, Pmed6, Pmed7, Pmed8]
Pstd_up = [Pmed1+Pstd1, Pmed2+Pstd2, Pmed3+Pstd3, Pmed4+Pstd4, Pmed5+Pstd5, Pmed6+Pstd6, Pmed7+Pstd7, Pmed8+Pstd8]
Pstd_low = [Pmed1-Pstd1, Pmed2-Pstd2, Pmed3-Pstd3, Pmed4-Pstd4, Pmed5-Pstd5, Pmed6-Pstd6, Pmed7-Pstd7, Pmed8-Pstd8]
Pmass_ind = [6.95, 7.45, 7.95, 8.45, 8.95, 9.45, 9.95, 10.1, 10.95]

two_Pmed = Pmed
two_Pstd_up = Pstd_up
two_Pstd_low = Pstd_low
two_Pmass_ind = Pmass_ind

readcol, 'data/PN_Mass_SFR/global_SFR_hopkins/global_SFR_mass_0.7_0.8.txt', P_mass1, P_SFR1, P_z1, format = 'D, D, D'

Pbit1 = (P_SFR1)[where(P_mass1 lt 7.5)]
Pbit2 = (P_SFR1)[where(P_mass1 ge 7.5 and P_mass1 lt 8.3)]
Pbit3 = (P_SFR1)[where(P_mass1 ge 8.3 and P_mass1 lt 8.8)]
Pbit4 = (P_SFR1)[where(P_mass1 ge 8.8 and P_mass1 lt 9.3)]
Pbit5 = (P_SFR1)[where(P_mass1 ge 9.3 and P_mass1 lt 9.8)]
Pbit6 = (P_SFR1)[where(P_mass1 ge 9.8 and P_mass1 lt 11.3)]

Pstd1 = stddev(Pbit1)
Pstd2 = stddev(Pbit2)
Pstd3 = stddev(Pbit3)
Pstd4 = stddev(Pbit4)
Pstd5 = stddev(Pbit5)
Pstd6 = stddev(Pbit6)

Pmed1 = median(Pbit1, /even)
Pmed2 = median(Pbit2, /even)
Pmed3 = median(Pbit3, /even)
Pmed4 = median(Pbit4, /even)
Pmed5 = median(Pbit5, /even)
Pmed6 = median(Pbit6, /even)

Pmed = [Pmed1, Pmed2, Pmed3, Pmed4, Pmed5, Pmed6]
Pstd_up = [Pmed1+Pstd1, Pmed2+Pstd2, Pmed3+Pstd3, Pmed4+Pstd4, Pmed5+Pstd5, Pmed6+Pstd6]
Pstd_low = [Pmed1-Pstd1, Pmed2-Pstd2, Pmed3-Pstd3, Pmed4-Pstd4, Pmed5-Pstd5, Pmed6-Pstd6]
Pmass_ind = [7.25, 7.9, 8.55, 9.05, 9.55, 10.55]

three_Pmed = Pmed
three_Pstd_up = Pstd_up
three_Pstd_low = Pstd_low
three_Pmass_ind = Pmass_ind

readcol, 'data/PN_Mass_SFR/global_SFR_hopkins/global_SFR_mass_0.8_1.1.txt', P_mass1, P_SFR1, P_z1, format = 'D, D, D'

Pbit1 = (P_SFR1)[where(P_mass1 lt 7.95)]
Pbit2 = (P_SFR1)[where(P_mass1 ge 7.95 and P_mass1 lt 8.45)]
Pbit3 = (P_SFR1)[where(P_mass1 ge 8.45 and P_mass1 lt 8.95)]
Pbit4 = (P_SFR1)[where(P_mass1 ge 8.95 and P_mass1 lt 9.45)]
Pbit5 = (P_SFR1)[where(P_mass1 ge 9.45 and P_mass1 lt 9.95)]
Pbit6 = (P_SFR1)[where(P_mass1 ge 9.95 and P_mass1 lt 11.95)]

Pstd1 = stddev(Pbit1)
Pstd2 = stddev(Pbit2)
Pstd3 = stddev(Pbit3)
Pstd4 = stddev(Pbit4)
Pstd5 = stddev(Pbit5)
Pstd6 = stddev(Pbit6)

Pmed1 = median(Pbit1, /even)
Pmed2 = median(Pbit2, /even)
Pmed3 = median(Pbit3, /even)
Pmed4 = median(Pbit4, /even)
Pmed5 = median(Pbit5, /even)
Pmed6 = median(Pbit6, /even)

Pmed = [Pmed1, Pmed2, Pmed3, Pmed4, Pmed5, Pmed6]
Pstd_up = [Pmed1+Pstd1, Pmed2+Pstd2, Pmed3+Pstd3, Pmed4+Pstd4, Pmed5+Pstd5, Pmed6+Pstd6]
Pstd_low = [Pmed1-Pstd1, Pmed2-Pstd2, Pmed3-Pstd3, Pmed4-Pstd4, Pmed5-Pstd5, Pmed6-Pstd6]
Pmass_ind = [7.7, 8.2, 8.7, 9.2, 9.7, 10.95]

four_Pmed = Pmed
four_Pstd_up = Pstd_up
four_Pstd_low = Pstd_low
four_Pmass_ind = Pmass_ind

first = create_struct('mass', one_Pmass_ind, 'SFR', one_Pmed, 'err_up', one_Pstd_up, 'err_low', one_Pstd_low)
sec = create_struct('mass', two_Pmass_ind, 'SFR', two_Pmed, 'err_up', two_Pstd_up, 'err_low', two_Pstd_low)
third = create_struct('mass', three_Pmass_ind, 'SFR', three_Pmed, 'err_up', three_Pstd_up, 'err_low', three_Pstd_low)
fourth = create_struct('mass', four_Pmass_ind, 'SFR', four_Pmed, 'err_up', four_Pstd_up, 'err_low', four_Pstd_low)

return, create_struct('first', first, 'sec', sec, 'third', third, 'fourth', fourth)

end
