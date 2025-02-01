# Debye-type relaxation

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexComp(E, T, f):
    numer = E**(-1)
    denom = 1 + T*f*(2j/2)
    comComp = numer/denom
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def complexComp(E, T, freq):
    strComp = []
    losComp = []
    for f in freq:
        s, l = calc_complexComp(E, T, f)
        strComp.append(s)
        losComp.append(l)
    complexComp = [strComp, losComp]
    return complexComp

def reqParams():
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    retardTime = viscosity/infMod
    return infMod, retardTime

def freqAxis(retardTime):
    centerFreq = 1 / retardTime
    freq = np.linspace(0, centerFreq*2, 10**2)
    scaledFreq = [f*retardTime for f in freq]
    freqAxis = [freq, scaledFreq]
    return freqAxis

if __name__=='__main__':
    # calcul1ating dynamic compliance and loss tangent
    infMod, retardTime = reqParams()
    param_text = '(Einf = {0:.1f} MPa, tau = {1:.1f} ms)'.format(infMod/10**6, retardTime*10**3)
    freqAxis = freqAxis(retardTime)

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = 'freq*tau'
    strComp = complexComp(infMod, retardTime, freq)[0]
    losComp = complexComp(infMod, retardTime, freq)[1]

    y1 = strComp
    y2 = losComp
    y_label = "J', "+'J" /Pa'
    ylim = [np.max(y1)*(-0.2), np.max(y1)*1.2]
    label1 = 'Storage compliance'
    label2 = 'Loss compliance'
    c1 = 'r'
    c2 = 'b'
    legend_loc='upper right'
    savefile = './png/Debye_relaxation.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Debye-type relaxation '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(freq, y1, c=c1, label=label1)
    ax.scatter(freq, y2, c=c2, label=label2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()