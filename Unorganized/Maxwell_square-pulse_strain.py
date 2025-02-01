# response of Maxwell model to square-pulse strain

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        modulus = float(input('Enter modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        modulus = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 5 kPa s): '))*10**3
    except ValueError:
        viscosity = 5*10**4
    try:
        p_width = float(input('Enter pulse width (ms) (default = 20 ms): '))*10**(-3)
    except ValueError:
        p_width = 20*10**(-3)
    relaxTime = viscosity/modulus
    return modulus, relaxTime, p_width

def func_pulseStrain(p_width):
    try:
        t_end = float(input('Enter t (pulse width < t) (ms) (default = 80 ms): '))*10**(-3)
    except ValueError:
        t_end = 80*10**(-3)
    num_t = 200
    duty_ratio = p_width/t_end
    tim = np.linspace(-t_end, t_end, num_t)
    tim0 = tim[:int(num_t*(1-duty_ratio)/2)]
    tim1 = tim[int(num_t*(1-duty_ratio)/2):int(num_t*(1+duty_ratio)/2)]
    tim2 = tim[int(num_t*(1+duty_ratio)/2):]
    strain0 = [0 for t in tim0]
    strain1 = [1 for t in tim1]
    strain2 = [0 for t in tim2]
    return t_end, tim0, tim1, tim2, strain0, strain1, strain2


def func_Maxwell1(modulus, relaxTime, tim, p_width):
    stress = [modulus*np.exp(-(t+p_width)/relaxTime) for t in tim]
    return stress

def func_Maxwell2(modulus, relaxTime, tau, tim):
    stress = [-modulus*(np.exp(tau/relaxTime)-np.exp(-tau/relaxTime))*np.exp(-t/relaxTime) for t in tim]
    return stress

if __name__=='__main__':
    E, T, tau = reqParams()
    te, tim0, tim1, tim2, strain0, strain1, strain2 = func_pulseStrain(tau)
    stress0 = [0 for t in range(len(tim0))]
    stress1 = func_Maxwell1(E, T, tim1, tau)
    stress2 = func_Maxwell2(E, T, tau, tim2)
    strain = strain0 + strain1 + strain2
    stress = stress0 + stress1 + stress2
    min_strain = np.min(strain)
    max_strain = np.max(strain)
    min_stress = np.min(stress)
    max_stress = np.max(stress)
    
    param_text = ' (E = {0:.1f} MPa, tau = {1:.1f} ms)'.format(E/10**6, T*10**3)

    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('Maxwell model for square pulse strain'+param_text)
    ax1.set_xlabel('time /s')
    ax1.set_ylabel('strain /')
    ax1.set_xlim(-te, te)
    ax1.set_ylim(np.max(strain)*(-0.2), np.max(strain)*1.2)
    ax1.grid()
    ax1.set_axisbelow(True)
    ax1.scatter(tim0, strain0, c='b')
    ax1.scatter(tim1, strain1, c='b')
    ax1.scatter(tim2, strain2, c='b')

    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('time /s')
    ax2.set_ylabel('stress /Pa')
    ax2.set_xlim(-te, te)
    ax2.set_ylim(abs(min_stress)*(-1.2), abs(max_stress)*1.2)
    ax2.grid()
    ax2.set_axisbelow(True)
    ax2.scatter(tim0, stress0, c='r')
    ax2.scatter(tim1, stress1, c='r')
    ax2.scatter(tim2, stress2, c='r')

    fig.savefig('./png/Maxwell_square-pulse_strain.png')

    plt.show()