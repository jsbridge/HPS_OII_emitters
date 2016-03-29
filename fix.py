#! /usr/bin/env python
import os
from glob import glob
import numpy as np

folders = glob('1*')+glob('2*')+glob('3*')+glob('4*')+glob('5*')+glob('6*')+glob('7*')+glob('8*')

for folder in folders:
    chi0, age0, gmass0, E0, smass0 = [], [], [], [], []
    chi1, age1, gmass1, E1, smass1 = [], [], [], [], []
    chi2, age2, gmass2, E2, smass2 = [], [], [], [], []
    chi3, age3, gmass3, E3, smass3 = [], [], [], [], []
    chi_med, age_med, gmass_med, E_med, smass_med = [], [], [], [], []
    os.chdir(folder)
    f = open('hps'+str(folder)+'_0.txt', 'r')
    for line in f:
    	col = line.split()
	chi0.append(float(col[1]))
	age0.append(float(col[2]))
	gmass0.append(float(col[3]))
	E0.append(float(col[4]))
	smass0.append(float(col[5]))
    chi_med.append(np.median(chi0))
    age_med.append(np.median(age0))
    gmass_med.append(np.median(gmass0))
    E_med.append(np.median(E0))
    smass_med.append(np.median(smass0))
    f = open('hps'+str(folder)+'_1.txt', 'r')
    for line in f:
    	col = line.split()
	chi1.append(float(col[1]))
	age1.append(float(col[2]))
	gmass1.append(float(col[3]))
	E1.append(float(col[4]))
	smass1.append(float(col[5]))
    chi_med.append(np.median(chi1))
    age_med.append(np.median(age1))
    gmass_med.append(np.median(gmass1))
    E_med.append(np.median(E1))
    smass_med.append(np.median(smass1))
    f = open('hps'+str(folder)+'_2.txt', 'r')
    for line in f:
    	col = line.split()
	chi2.append(float(col[1]))
	age2.append(float(col[2]))
	gmass2.append(float(col[3]))
	E2.append(float(col[4]))
	smass2.append(float(col[5]))
    chi_med.append(np.median(chi2))
    age_med.append(np.median(age2))
    gmass_med.append(np.median(gmass2))
    E_med.append(np.median(E2))
    smass_med.append(np.median(smass2))
    f = open('hps'+str(folder)+'_3.txt', 'r')
    for line in f:
    	col = line.split()
	chi3.append(float(col[1]))
	age3.append(float(col[2]))
	gmass3.append(float(col[3]))
	E3.append(float(col[4]))
	smass3.append(float(col[5]))
    chi_med.append(np.median(chi3))
    age_med.append(np.median(age3))
    gmass_med.append(np.median(gmass3))
    E_med.append(np.median(E3))
    smass_med.append(np.median(smass3))
    print folder, [('%.2f' % i) for i in chi_med], [('%.2f' % i) for i in age_med], [('%.2f' % i) for i in gmass_med], [('%.2f' % i) for i in E_med], [('%.2f' % i) for i in smass_med]
    os.chdir('..')
    
