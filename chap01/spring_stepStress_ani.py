# ordinary differential equation of spring (step stress) with animation
'''
バネ要素単独では常微分方程式は不要
'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

# variables
try:
    mod = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    mod = 2*10**5             # [Pa] modulus
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of spring width

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
el = e*l                        # [m] elongation

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.grid()
title_text = "spring (Hooke's elasticity): step stress, $\sigma_0$ = {0:.1f} MPa".format(s_i/10**6)
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)
x_rod1 = [0, l/4]
y = [0, 0]
ax.plot(x_rod1,y, c='b')

var_text = r'$\sigma_0$ = {0:.1f} MPa, $E$ = {1:.1f} MPa'.format(s_i/10**6,mod/10**6)
ax.text(0.6, 0.9, var_text, transform=ax.transAxes)
eq_text = r'$\epsilon$ = $\sigma_0/E$'
ax.text(0.6, 0.8, eq_text, transform=ax.transAxes)

rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, triangle, point, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod2 = [3*l/4 + el[i], l + el[i]]
    x_point = [0, l + el[i]]
    rod.set_data(x_rod2,y)
    x_tri = np.linspace(l/4, 3*l/4 + el[i],100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(el[i]+l/2)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    point.set_data(x_point,y)
    time_text.set_text(time_template % (t[i]))
    return rod, triangle, point, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/spring_stepStress_(mod={0:.1f}M,sigma={1:.1f}M).gif".format(mod/10**6,s_i/10**6)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()