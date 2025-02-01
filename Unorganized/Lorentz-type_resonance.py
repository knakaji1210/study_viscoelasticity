# Lorentz-type resonance

import numpy as np
import matplotlib.pyplot as plt

def calc_complexComp(E, T, f, f0):
    numer = (2j/2)*E**(-1)
    denom1 = 1 + T*(f-f0)*(2j/2)
    denom2 = 1 + T*(f+f0)*(2j/2)
    comComp = -numer/denom1 + numer/denom2
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def complexComp(E, T, freq, f0):
    strComp = []
    losComp = []
    for f in freq:
        s, l = calc_complexComp(E, T, f, f0)
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
    try:
        resoFreq = float(input('Enter resonant frequency (kHz) (default = 0.1 kHz): '))*10**3
    except ValueError:
        resoFreq = 0.1*10**3
    retardTime = viscosity/infMod
    return infMod, retardTime, resoFreq

def freqAxis(resoFreq, retardTime):
    freq = np.linspace(0, resoFreq*2, 10**2)
    scaledFreq = [f*retardTime for f in freq]
    freqAxis = [freq, scaledFreq]
    return freqAxis

if __name__=='__main__':
    # calcul1ating dynamic compliance
    infMod, retardTime, resoFreq = reqParams()
    param_text = '(Einf = {0:.1f} MPa, tau = {1:.1f} ms, F0 = {2:.1f} kHz)'.format(infMod/10**6, retardTime*10**3, resoFreq/10**3)
    freqAxis = freqAxis(resoFreq, retardTime)

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = 'freq*tau'
    strComp = complexComp(infMod, retardTime, freq, resoFreq)[0]
    losComp = complexComp(infMod, retardTime, freq, resoFreq)[1]

    y1 = strComp
    y2 = losComp
    y_label = "J', "+'J" /Pa'
    ylim = [np.min(y1)*1.5, np.max(y2)*1.2]
    label1 = 'Storage compliance'
    label2 = 'Loss compliance'
    c1 = 'r'
    c2 = 'b'
    legend_loc='upper right'
    savefile = './png/Lorentz_resonance.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Lorentz-type resonance '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(scaledFreq, y1, c=c1, label=label1)
    ax.scatter(scaledFreq, y2, c=c2, label=label2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()