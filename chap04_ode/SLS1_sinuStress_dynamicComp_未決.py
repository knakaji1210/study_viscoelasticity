# ordinary differential equation of SLS1 model (sinusoidal stress)
# freqSweep version -> storage & loss compliance

'''
How to use
% python3 SLS1_sinuStress_dynamicComp.py args[1]
args: -log
"-log"をつけると縦軸をログスケールに変換
何もついていないか、間違えたものがついている時はリニアスケールで表示
'''

import sys
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

if len(sys.argv) == 1:
    axisoption = ""
else:
    axisoption = sys.argv[1]
if axisoption == "-log":
    pass
else:
    axisoption = ""

'''
テキストの式(4.11)をベースに組み立てる
'''

def SLS1_sinuStress(e, t, samp, af, insMod, k, tau):
# e: strain, s: stress, infMod: relaxation modulus, tau: relaxation time
# ここではsampとafを指定し、この中でsの関数を作り振動歪みを実現
    s = samp*np.sin(af*t)
    dsdt = samp*af*np.cos(af*t)
    dedt = (k*s/insMod + tau*dsdt/insMod - e)/tau   # (4.11'')
    return dedt

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

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

# external sinusoidal stress
try:
    log_af_min = float(input('log(af_min) for forced oscillation (default=-2.0): '))
except ValueError:
    log_af_min = -2.0
try:
    log_af_max = float(input('log(af_max) for forced oscillation (default=1.0): '))
except ValueError:
    log_af_max = 1.0
try:
    num_freq = int(input('number of anglar frequency (default=30): '))
except ValueError:
    num_freq = 30
try:
    samp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))*10**6
except ValueError:
    samp = 0.2*10**6

af_list = np.logspace(log_af_min, log_af_max, num_freq)
aft_list = af_list*tau

dt = 0.01                  # [s] interval time
e0 = 0                     # [MPa] initial stress

s = []              # 周波数掃引の全ての入力信号（応力）を格納
e = []              # 周波数掃引の全ての出力信号（歪み）を格納
eamp = []           # 各周波数での出力振幅の最大値を格納
epha = []           # 各周波数での出力信号の位相を格納

for af in af_list:
    tmax = 10*np.pi/af          # [s] duration time
    dt = 0.01                   # [s] interval time
    t = np.arange(0, tmax, dt)  # time during sinusoidal stress
    s_af = samp*np.sin(af*t)    # 入力信号
    s_af_last = s_af[int(0.8*len(s_af)):]
    s.extend(s_af)            
    # solution of ODE
    sol = odeint(SLS1_sinuStress, e0, t, args=(samp,af,insMod,k,tau))
    e_af = sol[:,0]        # [] strain
    e_af_last = e_af[int(0.8*len(e_af)):]
    e_max = np.max(e_af_last)
    eamp.append(e_max)        
    ind = getNearestIndex2value(e_af_last,0)    # 出力信号が0になるindexを抽出
    e_pha = (180/np.pi)*np.arcsin(np.abs(s_af_last[ind])/samp)
    epha.append(e_pha)  

s = np.array(s)
e = np.array(e)

# J', J"の計算
'''
J' = epsilon*cos(theta)/samp
J" = epsilon*sin(theta)/same
'''
cos_epha = [np.cos(np.radians(epha)) for epha in epha]
sin_epha = [np.sin(np.radians(epha)) for epha in epha]
if axisoption == "-log":    
    strComp = [e*cosp/samp for (e,cosp) in zip(eamp,cos_epha)]
    losComp = [e*sinp/samp for (e,sinp) in zip(eamp,sin_epha)]
else:
    strComp = [e*cosp/samp for (e,cosp) in zip(eamp,cos_epha)]
    losComp = [e*sinp/samp for (e,sinp) in zip(eamp,sin_epha)]

# scaling for figure
s = s/10**6                   # MPaスケール

eamp_max = np.max(eamp)

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax1.grid()
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
title_text = "SLS1 model: sinusoidal stress"
fig.suptitle(title_text)
ax1.set_axisbelow(True)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\omega_f$ $\tau$')
if axisoption == "-log":
    ax1.set_ylim(10**(-2), 10**(3))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /Pa$^{{{-1}}}$')
    ax2.set_ylim(10**(-2), 10**(3))
    ax2.set_ylabel(r'loss compliance, $J^{{\prime\prime}}$ /Pa$^{{{-1}}}$')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
else:
#    ax1.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /MPa$^{{{-1}}}$')
#    ax2.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax2.set_ylabel(r'loss complianca, $J^{{\prime\prime}}$ /MPa$^{{{-1}}}$')

var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(samp/10**6,E1/10**6,E2/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($k$$\sigma$/$E_i$ + $\tau$/$E_i$ d$\sigma$/d$t$ - $\epsilon$)/$\tau$'
ax1.text(0.1, 0.7, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax1.text(0.1, 0.6, res_text, transform=ax2.transAxes)

ax1.plot(aft_list,eamp, 'ro-', label=r'$J^{{\prime}}$')
ax1.set_xscale('log')
ax2.plot(aft_list,epha, 'bo-', label=r'$J^{{\prime\prime}}$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

if axisoption == "-log":
    savefile = "./png/SLS1_sinuStress_dynamicComp(log)_(simga={0:.2f}MPa,tau={1:.1f}s,E1={2:.1f}MPa,E2={3:.1f}MPa,eta={4:.1f}kPas).png".format(samp/10**6,tau,E1/10**6,E2/10**6,eta/10**3)
else:
    savefile = "./png/SLS1_sinuStress_dynamicComp(linear)_(sigma={0:.2f}MPa,tau={1:.1f}s,E1={2:.1f}MPa,E2={3:.1f}MPa,eta={4:.1f}kPas).png".format(samp/10**6,tau,E1/10**6,E2/10**6,eta/10**3)

fig.savefig(savefile, dpi=300)

plt.show()