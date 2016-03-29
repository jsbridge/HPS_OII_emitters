import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def power(x, a, beta):
    return a * x**beta

def line(x, b, m):
    return m*x + b


def size_mass():

    size, mass = np.loadtxt('morphology/OII_sizes.txt', unpack = True, usecols = (3, 4))

    size = np.log10(size)


    popt, pcov = curve_fit(line, mass, size)
    alpha = popt[1]
    perr = np.sqrt(np.diag(pcov))
    
    print alpha
    print perr[1]

    plt.plot(mass, size, 'k*')
    plt.plot(mass, line(mass, popt[0], popt[1]), 'b-')
    plt.show()

    return None
