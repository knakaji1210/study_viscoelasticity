# ordinary differential equation of Voigt model (step stress)  with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

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
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of dashpot width

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
el = e*l                                    # [m] elongation
dedt = np.array([0.0]+[(e[k+1]-e[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt
s_s = E*e                                   # stress on spring
s_d = eta*dedt                              # stress on dashpot

# scaling for figure
s = s/10**6                     # MPaスケール
s_s = s_s/10**6                 # MPaスケール
s_d = s_d/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Voigt model: step stress"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([0, l/4],y_0, c='b')
ax.plot([l/4,l/4],[1,-1], c='b')
ax.plot([l/4, l/2],y_1, c='b')
ax.plot([l/4, l/2],y_2, c='b')
ax.plot(0,0,'ro', markersize='10')
# for dashpot
x_d1 = [l/2, (0.92+1/4)*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([l/2,l/2],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=(l/2, -w+1), width=2*l/3, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\sigma_0$ = {0:.1f} MPa, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(s_i/10**6,E/10**6,eta/10**3)
ax.text(0.6, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma_0$/$E$ - $\epsilon$)/$\tau$'
ax.text(0.6, 0.8, eq_text, transform=ax.transAxes)
tau = eta/E
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.6, 0.7, res_text, transform=ax.transAxes)

# for common
var, = ax.plot([],[], 'b', animated=True)
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# for spring
rod_sp, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
# for dashpot
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)

# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return var, rod, point, rod_sp, triangle, rod_da, damper, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = l + el[i]
    x_c = l/4 + x_o
    var.set_data([x_c,x_c],[1,-1])
    x_rod = [x_c, l/4+x_c]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_c],[0])
    # for spring
    x_rod_sp = [x_o, x_c]
    rod_sp.set_data(x_rod_sp,y_2)
    x_tri = np.linspace(l/2, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/2)/(el[i]+l/2)-np.pi/2+0.1))-3)
    triangle.set_data(x_tri,y_tri)
    # for dashpot
    x_rod_da = [x_o-l/4-0.1*l, x_c]
    x_damp = x_o-l/4-0.1*l      # 0.1*lは位置調整のため
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (t[i]))
    return var, rod, point, rod_sp, triangle, rod_da, damper, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/Voigt_stepStress_(sigma={0:.1f}MPa,mod={1:.1f}MPa,eta={2:.1f}kPas,).gif".format(s_i/10**6,E/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()