# ordinary differential equation of SLS1 model (step strain)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(4.11)をベースに組み立てる
'''

def SLS1_stepStrain(s, t, e, infMod, k, tau):
# e: strain, s: stress, infMod: relaxation modulus, tau: relaxation time
# ここでは下でargsとしてe=e0を入れてステップ応力を実現
    dsdt = k*(infMod*e - s)/tau   # (4.11)
    return dsdt

# variables
try:
    E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
except ValueError:
    E1 = 10**6                    # [Pa] modulus
try:
    E2 = float(input('modulus 2 [MPa] (default = 0.5 MPa): '))*10**6
except ValueError:
    E2 = 5*10**5                  # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 1000.0 kPa s): '))*10**3
except ValueError:
    eta = 10**6                   # [Pa s] viscosity
tau = eta/E2                  # [s] retardation time
insMod = E1                   # [Pa] instantaneous modulus
infMod = E1*E2/(E1+E2)        # [Pa] equilibrium modulus
k = insMod/infMod

# initial condition
try:
    e0 = float(input('step strain [] (default = 0.2): '))
except ValueError:
    e0 = 0.2                # [] step strain
s0 = E1*e0                  # [Pa] initial stress

tmax = 10                   # [s] duration time
dt = 0.01                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
e = np.concatenate([zeros,ones*e0])

# solution of ODE
sol = odeint(SLS1_stepStrain, s0, t_a, args=(e0,infMod,k,tau))
s = np.concatenate([zeros,sol[:,0]])        # [] strain
e1 = s/E1                                   # [] strain of spring1
e2 = e - e1                                 # [] strain of Voigt component
s1 = E2*e2                                  # [MPa] stress on sping in Voigt component
de2dt = np.array([0.0]+[(e2[k+1]-e2[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt (Voigt部分)
s2 = eta*de2dt                              # [MPa] stress on dashpot
sv = s1 + s2                                # [MPa] Voigt部分の応力確認用

# scaling for figure
s = s/10**6                   # MPaスケール
s1 = s1/10**6                 # MPaスケール
s2 = s2/10**6                 # MPaスケール
sv = sv/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "SLS1 model: step strain"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$ /MPa')
ax1.set_ylim(-0.1*np.max(e),1.5*np.max(e))
ax2.set_ylim(-0.1*np.max(s),1.5*np.max(s))

var_text = r'$\epsilon_0$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(e0,E1/10**6,E2/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = ($k$/$\tau$)($E_{{{\infty}}}$$\epsilon_0$ - $\sigma$)'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s, $\tau$/$k$ = {1:.1f} s'.format(tau, tau/k)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

ax1.plot(t, e, 'b', label='$\epsilon$ (input)')
ax1.plot(t, e1, 'g', ls="dashed", label='$\epsilon$ (spring1)')
ax1.plot(t, e2, 'y', ls="dashed", label='$\epsilon$ (Voigt)')
ax1.legend(loc='upper right')
ax2.plot(t, s, 'r', label='$\sigma$ (output)')
ax2.plot(t, s1, 'g', ls="dashed", label='$\sigma$ (spring2)')
ax2.plot(t, s2, 'y', ls="dashed", label='$\sigma$ (dashpot)')
#ax2.plot(t, sv, 'k', ls="dashed", label='$\sigma$ (Voigt)')
ax2.legend(loc='upper right')


savefile = "./png/SLS1_stepStrain_(sigma={0:.2f}MPa,E1={1:.1f}MPa,E2={2:.1f}MPa,eta={3:.1f}kPas).png".format(s0/10**6,E1/10**6,E2/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()