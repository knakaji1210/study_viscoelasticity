# ordinary differential equation of dashpot (step stress) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

def dashpot(e, t, s, eta):
    dedt = s/eta        # e: strain, s: step stress, eta: viscosity
    return dedt

# variables
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5           # [Pa s] viscosity
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of dashpot width

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

# solution of ODE
sol = odeint(dashpot, e0, t_a, args=(s_i,eta))
e = np.concatenate([zeros,sol[:,0]])            # [] strain
el = e*l                                        # [m] elongation

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.grid()
title_text = "dashpot (Newton's viscosity): step stress, $\sigma_0$ = {0:.1f} MPa".format(s_i/10**6)
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)
x_rod1 = [0, l/4]
y = [0, 0]
x_d1 = [l/4, 0.92*l]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_rod1,y, c='b')
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([l/4,l/4],[w,-w], c='b')
rect = patches.Rectangle(xy=(l/4, -w), width=2*l/3, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\sigma_0$ = {0:.1f} MPa, $\eta$ = {1:.1f} kPa s'.format(s_i/10**6,eta/10**3)
ax.text(0.6, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma_0$/$\eta$'
ax.text(0.6, 0.8, eq_text, transform=ax.transAxes)

rod, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, damper, point, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod2 = [l/2 + el[i], l + el[i]]
    x_damp = (l/2)+el[i]
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    x_point = [0, l + el[i]]
    rod.set_data(x_rod2,y)
    damper.set_data(x_damper,y_damper)
    point.set_data(x_point,y)
    time_text.set_text(time_template % (t[i]))
    return rod, damper, point, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/dashpot_stepStress_ani_(eta={0:.1f}k,sigma={1:.1f}M).gif".format(eta/10**3,s_i/10**6)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()