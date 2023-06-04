# ordinary differential equation of Voigt model (step stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(1.22)をベースに組み立てる
'''

def Voigt(e, t, s_i, E, eta):
    tau = eta/E                 # retardation time [s]
    e_inf = s_i/E
    dedt = (e_inf - e)/tau      # (1.22)
    return dedt

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
    s_i = float(input('step stress [MPa] (default = 0.1 MPa): '))*10**6
except ValueError:
    s_i = 10**5             # [Pa] stepp stress
e0 = 0                      # [] initial strain

tmax = 10                   # [s] duration time
dt = 0.05                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
s = np.concatenate([zeros,ones*s_i])

# solution of ODE
sol = odeint(Voigt, e0, t_a, args=(s_i,E,eta))
e = np.concatenate([zeros,sol[:,0]])        # [] strain
dedt = np.array([0.0]+[(e[k+1]-e[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt
s_s = E*e                                   # stress on spring
s_d = eta*dedt                              # stress on dashpot

# scaling for figure
s = s/10**6                     # MPaスケール
s_s = s_s/10**6                 # MPaスケール
s_d = s_d/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(111, xlabel='$t$ /s')
ax2 = ax1.twinx()
ax1.grid()
ax2.grid(ls='dotted')
title_text = "Voigt model: step stress"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$')
ax1.set_ylim(-0.1*np.max(e),1.5*np.max(e))
ax2.set_ylim(-0.1*np.max(s),1.5*np.max(s))

var_text = r'$\sigma_0$ = {0:.1f} MPa, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(s_i/10**6,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma_0$/$E$ - $\epsilon$)/$\tau$'
ax1.text(0.1, 0.8, eq_text, transform=ax1.transAxes)
tau = eta/E
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax1.text(0.1, 0.7, res_text, transform=ax1.transAxes)

ax1.plot(t, e, 'b', label='$\epsilon$ (output)')
ax1.legend(loc='upper right')
ax2.plot(t, s, 'r', label='$\sigma$ (input)')
ax2.plot(t, s_s, 'g', ls="dashed", label='$\sigma$ (spring)')
ax2.plot(t, s_d, 'y', ls="dashed", label='$\sigma$ (dashpot)')
ax2.legend(loc='lower right')

savefile = "./png/Voigt_stepStress_(sigma={0:.1f}MPa,mod={1:.1f}MPa,eta={2:.1f}kPas,).png".format(s_i/10**6,E/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()