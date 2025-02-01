# transient response of Voigt model

import numpy as np
import matplotlib.pyplot as plt

def sinusoidalExcitation(f, amp, time): # f must be angular frequency
    inputReal = []
    inputImag = []
    for t in time:
        inputSignal = np.exp((2j/2)*f*t)*amp
        inR = inputSignal.real
        inI = inputSignal.imag
        inputReal.append(inR)
        inputImag.append(inI)
    return inputReal, inputImag

def complexComp(E, T, f): # f must be angular frequency
    numer = E**(-1)
    denom = 1 + T*f*(2j/2)
    comComp = numer/denom
    return comComp

def transientResponse(E, T, f, amp, time):
    outputReal = []
    outputImag = []
    for t in time:
        outputSignal = complexComp(E, T, f)*(np.exp((2j/2)*f*t)-np.exp(-t/T))*np.heaviside(t,0.5)*amp
        outR = outputSignal.real
        outI = outputSignal.imag
        outputReal.append(outR)
        outputImag.append(outI)
    return outputReal, outputImag

def reqParams():
    try:
        freq = float(input('Enter vibration freq (Hz) (default = 1 Hz): '))
    except ValueError:
        freq = 1
    try:
        vibAmp = float(input('Enter vibration amplitude (MPa) (default = 1.0 MPa): '))*10**6
    except ValueError:
        vibAmp = 10**6
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        infMod = 10**5
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    retardTime = viscosity/infMod
    return freq, vibAmp, infMod, retardTime

def timeAxis(T):
    time = np.linspace(0, T*10, 4096)
    scaledTime = [t/T for t in time]
    timeAxis = [time, scaledTime]
    return timeAxis

if __name__=='__main__':
    # calcul1ating transientt response
    freq, vibAmp, infMod, retardTime = reqParams()
    angFreq = 2*np.pi*freq
    param_text = '(Einf = {0:.1f} MPa, tau = {1:.1f} ms, f = {2:.1f} Hz)'.format(infMod/10**6, retardTime*10**3, freq)
    time = timeAxis(retardTime)[0]
    scaledTime = timeAxis(retardTime)[1]
    time_label = 'time/tau'
    inputReal, inputImag = sinusoidalExcitation(angFreq, vibAmp, time)
    outputReal, outputImag = transientResponse(infMod, retardTime, angFreq, vibAmp, time)
    try:
        select = int(input('Selection (real part: 0, imaginary part: 1): '))
    except ValueError:
        select = 0

    # drawing graphs
    if select == 0:
        y1 = inputReal
        y2 = outputReal
        y1_label = 'stress /MPa'
        y2_label = 'strain /'
        y1lim = [np.min(y1)*3.0, np.max(y1)*3.0]
        y2lim = [np.min(y2)*1.5, np.max(y2)*1.5]
        label1 = 'input signal (real part = cos)'
        label2 = 'output signal'
        c1 = 'r'
        c2 = 'b'
        legend_loc='upper right'
        savefile = './png/Voigt_transient_response_real.png'

    if select == 1:
        y1 = inputImag
        y2 = outputImag
        y1_label = 'stress /MPa'
        y2_label = 'strain /'
        y1lim = [np.min(y1)*3.0, np.max(y1)*3.0]
        y2lim = [np.min(y2)*1.5, np.max(y2)*1.5]
        label1 = 'input signal (imaginary part = sin)'
        label2 = 'output signal'
        c1 = 'r'
        c2 = 'b'
        legend_loc='upper right'
        savefile = './png/Voigt_transient_response_imag.png'

    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax1 = fig.add_subplot(111)
    ax1.set_title('Voigt model '+param_text)
    ax1.set_ylim(y1lim[0], y1lim[1])
    ax1.set_xlabel('t/tau')
    ax1.set_ylabel(y1_label)
    ax1.plot(scaledTime, y1, c=c1, label=label1)

    ax2 = ax1.twinx()
    ax2.plot(scaledTime, y2, c=c2, label=label2)
    ax2.set_ylim(y2lim[0], y2lim[1])
    ax2.set_ylabel(y2_label)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc='upper right')
    ax1.grid()
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)    
    
    fig.savefig(savefile)  

    plt.show()