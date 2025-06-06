import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

infile = open("deltaW.out", "r")

qa = []
qc = []
pc = []
mu = []
ep = []

g1 = []
g2 = []
w1 = []
w2 = []
f1 = []
f2 = []

next(infile)

for line in infile: 

    numbers = line.split()

    c1 = float(numbers[0])
    c2 = float(numbers[1])
    c3 = float(numbers[2])
    c4 = float(numbers[3])
    c5 = float(numbers[4])

    qc.append(c1)
    qa.append(c2)
    pc.append(c3)
    mu.append(c4)
    ep.append(c5)

    c6  = float(numbers[8])
    c7  = float(numbers[9])
    c8  = float(numbers[10])
    c9  = float(numbers[11])
    c10 = float(numbers[12])
    c11 = float(numbers[13])

    g1.append(c6)
    w1.append(c7)
    f1.append(c8)
    g2.append(c9)
    w2.append(c10)
    f2.append(c11)

b0 = []
for p in pc:
    b0.append(0.08*p)

pmin = 0.
pmax = 1.4e-2
ymax = 1.e-6
    
fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'Circular: q_0=0.8: q_a=2.75, epsilon_a=0.2')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

ax = plt.subplot (3, 1, 1)

plt.xlim (pmin, pmax)

plt.plot    (b0, g1, color = 'blue', linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\gamma_1$')
plt.plot    (b0, g2, color = 'red',  linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\gamma_2$')

plt.axvline (0.1674*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,            color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$\gamma_1, \gamma_2$ (kHz)", fontsize = "15")
plt.legend (fontsize = '15')

inset_ax = inset_axes(ax, width="50%", height="30%", bbox_to_anchor=(0.1, 0., 1., 1.), bbox_transform=ax.transAxes, loc='upper center')

inset_ax.tick_params(axis='both', labelsize=10)

inset_ax.set_xlim (0.00, 0.014)
inset_ax.set_ylim (-0.5, 120.)

inset_ax.plot    (b0, g1, color = 'blue', linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
inset_ax.plot    (b0, g2, color = 'red',  linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)

inset_ax.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.subplot (3, 1, 2)

plt.xlim (pmin, pmax)

plt.plot    (b0, w1, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\omega_1$')
plt.plot    (b0, w2, color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\omega_2$')

plt.axvline (0.1674*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,            color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$\omega_1, \omega_2$ (kHz)", fontsize = "15")
plt.legend (fontsize = '15')

plt.subplot (3, 1, 3)

plt.xlim (pmin, pmax)

plt.plot    (b0, f1, color = 'blue', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$f_1$')
plt.plot    (b0, f2, color = 'red',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$f_2$')

plt.axvline (0.1674*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,            color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$f_1, f_2$", fontsize = "15")
plt.legend (fontsize = '15')

plt.tight_layout ()

#plt.show ()    
plt.savefig("gamma.pdf")
