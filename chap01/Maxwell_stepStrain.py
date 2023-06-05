# ordinary differential equation of Maxwell model (step strain)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(1.12)をベースに組み立てる
'''

def Maxwell(s, t, E, eta):
    tau = eta/E                 # retardation time [s]
    dsdt = -s/tau               # (1.12)
    return dsdt

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] viscosity

# initial condition
try:
    e_i = float(input('step strain [] (default = 0.4): '))
except ValueError:
    e_i = 0.4               # [] step strain
s0 = E*e_i                  # [Pa] initial stress

tmax = 10                   # [s] duration time
dt = 0.01                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
e = np.concatenate([zeros,ones*e_i])

# solution of ODE
sol = odeint(Maxwell, s0, t_a, args=(E,eta))
s = np.concatenate([zeros,sol[:,0]])        # [] stress
integral_s = np.array([s[:k+1].sum()*dt for k in range(len(s))])     # 簡易的なsの積分
e_s = s/E                                   # strain on spring
e_d = integral_s/eta                              # strain on dashpot

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,8))
title_text = "Maxwell model: step strain"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$')
ax1.set_ylim(-0.1*np.max(e),1.5*np.max(e))
ax2.set_ylim(-0.1*np.max(s),1.5*np.max(s))

var_text = r'$\epsilon_0$ = {0:.1f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(e_i,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
tau = eta/E
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

ax1.plot(t, e, 'b', label='$\epsilon$ (input)')
ax1.plot(t, e_s, 'g', ls="dashed", label='$\epsilon$ (spring)')
ax1.plot(t, e_d, 'y', ls="dashed", label='$\epsilon$ (dashpot)')
ax1.legend(loc='upper right')
ax2.plot(t, s, 'r', label='$\sigma$ (output)')
ax2.legend(loc='upper right')

savefile = "./png/Maxwell_stepStrain_(epsilon={0:.2f},mod={1:.1f}MPa,eta={2:.1f}kPas).png".format(e_i,E/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()