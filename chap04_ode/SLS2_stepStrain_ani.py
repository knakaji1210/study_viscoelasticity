# ordinary differential equation of SLS2 model (step strain) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

'''
テキストの式(4.37)をベースに組み立てる
'''

def SLS2_stepStrain(s, t, e, infMod, tau):
# e: strain, s: stress, infMod: relaxation modulus, tau: relaxation time
# ここでは下でargsとしてe=e0を入れてステップ応力を実現
    dsdt = (infMod*e - s)/tau   # (4.11)
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
tau = eta/E1                   # [s] retardation time
insMod = E1+E2                 # [Pa] instantaneous modulus
infMod = E2                    # [Pa] equilibrium modulus
k = insMod/infMod
l = 0.1                         # [m] equilibrium length
w = 0.5                         # ratio of dashpot width

# initial condition
try:
    e0 = float(input('step strain [] (default = 0.2): '))
except ValueError:
    e0 = 0.2                # [] step strain
s0 = insMod*e0              # [Pa] initial stress

tmax = 10                   # [s] duration time
dt = 0.05                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
e = np.concatenate([zeros,ones*e0])
# 直列なのでl0=2lとなっている
el = e*2*l                                  # [m] elongation

# solution of ODE
sol = odeint(SLS2_stepStrain, s0, t_a, args=(e0,infMod,tau))
s = np.concatenate([zeros,sol[:,0]])        # [] strain
s2 = E2*e                                   # [MPa] stress of spring2
s1 = s - s2                                 # [MPa] stress of Maxwell component
e1 = s1/E1                                  # [] strain of spring1
e2 = e - e1                                 # [] strain of dashpot
de2dt = np.array([0.0]+[(e2[k+1]-e2[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt (dashpot部分)
s1_c = eta*de2dt                            # [MPa] stress on dashpot（確認用）
el_s = e1*2*l                              # [m] elongation of spring1
el_d = e2*2*l                              # [m] elongation of dashpot

# scaling for figure
s = s/10**6                   # MPaスケール
s1 = s1/10**6                 # MPaスケール
s2 = s2/10**6                 # MPaスケール
s1_c = s1_c/10**6             # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "SLS2 model: step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.4)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([l/4, l/4+0.08*l],y_1, c='b')
ax.plot(0,0,'ro', markersize='10')
ax.plot([0, l/4],y_0, c='b')
ax.plot([l/4,l/4],[1,-1], c='b')
ax.plot([l/4, l],y_2, c='b')            # 7*l/8ではないかも
ax.plot([l/4,9*l/4],[-3,-3], c='g')
ax.plot([l/4,l/4],[-2.8,-3.2], c='g')
ax.plot([9*l/4,9*l/4],[-2.8,-3.2], c='g')
# for dashpot
x_d1 = [l/4+0.08*l, l/4+0.92*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([l/4+0.08*l,l/4+0.08*l],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=(l/4+0.08*l, -w+1), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\epsilon_0$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(e0,E1/10**6,E2/10**6,eta/10**3)
ax.text(0.5, 0.9, var_text, transform=ax.transAxes)
#eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$'
#ax.text(0.5, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.5, 0.7, res_text, transform=ax.transAxes)
ax.text(0.4, 0.15, '$l_0$', transform=ax.transAxes)

# for common
bar, = ax.plot([],[], 'b', animated=True)
rod1, = ax.plot([],[], 'b', animated=True)
rod2, = ax.plot([],[], 'b', animated=True)
rod3, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# for spring1
triangle1, = ax.plot([],[], 'b', animated=True)
# for dashpot
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# for spring2
triangle2, = ax.plot([],[], 'b', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod1, rod2, rod3, point, triangle1, rod_da, damper, triangle2, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = 7*l/4 + el[i]
    bar.set_data([l/2+x_o,l/2+x_o],[1,-1])
    x_rod1 = [l/2+x_o, 3*l/4+x_o]
    x_rod2 = [l/4+x_o, l/2+x_o]
    x_rod3 = [-l/4+x_o, l/2+x_o]
    rod1.set_data(x_rod1,y_0)
    rod2.set_data(x_rod2,y_1)
    rod3.set_data(x_rod3,y_2)
    point.set_data([3*l/4+x_o],[0])
    # for dashpot
    x_da = 3*l/4 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (t[i]))
    # for spring1
    x_sp1 = x_o-el_s[i]-l/4
    x_tri1 = np.linspace(x_da+3*l/4, l/4+x_o,100)
    y_tri1 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri1 - x_da-3*l/4)/(x_o-x_da-l/2)-np.pi/2+0.1))+1)
    triangle1.set_data(x_tri1,y_tri1)
    # for spring2
    x_sp2 = -l/4+x_o
    x_tri2 = np.linspace(l, x_sp2,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - x_sp2)/(x_sp2-l)-np.pi/2+0.1))-3)
    triangle2.set_data(x_tri2,y_tri2)
    return bar, rod1, rod2, rod3, point, triangle1, rod_da, damper, triangle2, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/SLS2_stepStrain_(epsilon={0:.2f},E1={1:.1f}MPa,E2={2:.1f}MPa,eta={3:.1f}kPas).gif".format(e0,E1/10**6,E2/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()