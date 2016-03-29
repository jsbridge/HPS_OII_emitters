#! /usr/bin/env python
#
# This script reads all the outputs from morphcode for COSMOS and GOODS and plots Gini vs. M20
#
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

folders, flag, mag, G, M20, flag1, SN, R_circ = [], [], [], [], [], [], [], []

os.chdir('cosmos')
folders1 = glob('1*')+glob('2*')+glob('3*')

for i, folder in enumerate(folders1):
    os.chdir(folder)
    f = open('outmorph.dat', 'r')
    for line in f:
    	if "#" in line:
            continue
    	col = line.split()
	G.append(float(col[19]))
	M20.append(float(col[20]))
        flag.append(float(col[21]))
        SN.append(col[5])
        R_circ.append(col[6])
        if (col[21] == '1'):
            #print col[20], col[19], folder
            flag1.append(1)
    f.close()
    os.chdir('..')
    if float(col[5]) < 2.5:
        print 'S/N < 2.5:', folder, col[5]
    folders.append(folder)

    

os.chdir('../goods')
folders2 = glob('3*')+glob('4*')

for i, folder in enumerate(folders2):
    os.chdir(folder)
    f = open('outmorph.dat', 'r')
    for line in f:
    	if "#" in line:
            continue
    	col = line.split()
	G.append(float(col[19]))
	M20.append(float(col[20]))
        flag.append(float(col[21]))
        SN.append(col[5])
        R_circ.append(col[6])
        if col[21] == '1':
            #print col[20], col[19], folder
            flag1.append(1)
    f.close()
    os.chdir('..')
    if float(col[5]) < 2.5:
        print 'S/N < 2.5:', folder, col[5]
    folders.append(folder)

os.chdir('..')

# Write the stuff to a file, in case it's needed
f = open('G_M20_withflag.txt', 'w')
outstring = zip(folders, G, M20, R_circ, flag)
for line in outstring:
    f.write("    ".join(str(x) for x in line) + '\n')
f.close()
print len(G), len(M20)

# Take out the flag = 1 galaxies
b = np.array(flag) < 0.5
G = b*G
M20 = b*M20
G = np.delete(G, np.where(G == 0))
M20 = np.delete(M20, np.where(M20 == 0))

print 'Number of galaxies with flag = 1:',len(flag1)


noise1 = np.random.normal(1, 0.003, len(G))
G = G * noise1
noise2 = np.random.normal(1, 0.003, len(G))
M20 = M20 * noise2

minorLocator   = MultipleLocator(5)

fig = plt.figure()
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.4)

plt.scatter(M20, G, color = 'green', linewidth = 1.4)
ax.xaxis.set_minor_locator(minorLocator)
plt.plot([-0.5, -2.75], [0.43, 0.655], 'k--', linewidth = 1.4)
plt.vlines(-0.35, 0.685, 0.615, linewidth = 1.4)
plt.hlines(0.65, -0.25, -0.45, linewidth = 1.4)
plt.xlabel('M$_{20}$')
plt.ylabel('G')
plt.ylim(0.37, 0.7)
plt.xlim(-3,0)
#plt.gca().xaxis.set_major_locator(MaxNLocator(prune='upper'))
plt.gca().invert_xaxis()
plt.xticks([0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0], ['0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0'])
plt.gca().tick_params(width = 1.3)
plt.savefig('GinivsM20.pdf')
plt.close()
plt.show()

# Write the stuff to a file, in case it's needed
f = open('G_M20_noflag.txt', 'w')
outstring = zip(folders, G, M20, R_circ, size)
for line in outstring:
    f.write("    ".join(str(x) for x in line) + '\n')
f.close()
