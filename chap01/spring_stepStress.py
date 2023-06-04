# ordinary differential equation of spring (step stress)
'''
バネ要素単独では常微分方程式は不要
'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# variables
try:
    mod = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    mod = 2*10**5             # [Pa] modulus

# initial condition
try:
    s_i = float(input('step stress [MPa] (default = 0.1 MPa): '))*10**6
except ValueError:
    s_i = 10**5             # [Pa] stepp stress
e0 = 0                      # [] initial strain

tmax = 2                    # [s] duration time
dt = 0.05                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.5*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
s = np.concatenate([zeros,ones*s_i])

# solution of ODE（ここではそれをする必要はない）
e = s/mod

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(111, xlabel='$t$ /s')
ax2 = ax1.twinx()
ax1.grid()
ax2.grid(ls='dotted')
title_text = "spring (Hooke's elasticity): step stress, $\sigma_0$ = {0:.1f} MPa".format(s_i/10**6)
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$')
ax1.set_ylim(-0.1*np.max(e),1.5*np.max(e))
ax2.set_ylim(-0.1*np.max(s),1.5*np.max(s))

var_text = r'$\sigma_0$ = {0:.1f} MPa, $E$ = {1:.1f} MPa'.format(s_i/10**6,mod/10**6)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'$\epsilon$ = $\sigma_0/E$'
ax1.text(0.1, 0.8, eq_text, transform=ax1.transAxes)

ax1.plot(t, e, 'b', label='$\epsilon$ (output)')
ax1.legend(loc='upper right')
ax2.plot(t, s, 'r', label='$\sigma$ (input)')
ax2.legend(loc='lower right')

savefile = "./png/spring_stepStress_(mod={0:.1f}M,sigma={1:.1f}M).png".format(mod/10**6,s_i/10**6)
fig.savefig(savefile, dpi=300)

plt.show()