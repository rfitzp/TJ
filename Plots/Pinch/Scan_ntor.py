# Scan_ntor.py

# Plots growth rates versus toroidal mode number

import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("../../Outputs/Pinch/mm1.out", "r")

nn = []
g1 = []
g2 = []
g3 = []

for line in infile: 

    numbers = line.split()
    q0      = float(numbers[2])
    qa      = float(numbers[3])
    c1      = float(numbers[4])
    c2      = float(numbers[5])
    c3      = float(numbers[6])
    c4      = float(numbers[7])
    nn.append(c1)
    g1.append(c2)
    g2.append(c3)
    g3.append(c4)

infile = open("../../Outputs/Pinch/m0.out", "r")

nn0 = []
g10 = []
g20 = []
g30 = []

for line in infile: 

    numbers = line.split()
    c1      = float(numbers[4])
    c2      = float(numbers[5])
    c3      = float(numbers[6])
    c4      = float(numbers[7])
    nn0.append(c1)
    g10.append(c2)
    g20.append(c3)
    g30.append(c4)

infile = open("../../Outputs/Pinch/m1.out", "r")    
    
nn1 = []
g11 = []
g21 = []
g31 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[4])
    c2      = float(numbers[5])
    c3      = float(numbers[6])
    c4      = float(numbers[7])
    nn1.append(c1)
    g11.append(c2)
    g21.append(c3)
    g31.append(c4)    

fig = plt.figure (figsize = (8., 8.))
fig.canvas.manager.set_window_title (r'Pinch Code: Toroidal Mode Number Scan')
plt.rc ('xtick', labelsize=15) 
plt.rc ('ytick', labelsize=15)

plt.subplot (3, 1, 1)

plt.plot (nn, g1, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'o', markersize = 5, label = r"$d/a=0.01$")
plt.plot (nn, g2, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'v', markersize = 5, label = r"$d/a=0.05$")
plt.plot (nn, g3, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 's', markersize = 5, label = r"$d/a=0.10$")

plt.axvline (-qa, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.title  ("$m=-1$",            fontsize = "15")
plt.xlabel (r'$n\,\epsilon_a$',  fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (3, 1, 2)

plt.plot (nn0, g10, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'o', markersize = 5, label = r"$d/a=0.01$")
plt.plot (nn0, g20, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'v', markersize = 5, label = r"$d/a=0.05$")
plt.plot (nn0, g30, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 's', markersize = 5, label = r"$d/a=0.10$")

plt.axhline (0., color = 'black', linewidth=1.5, linestyle = 'dotted')

plt.title  ("$m=0$", fontsize = "15")
plt.xlabel (r'$n\,\epsilon_a$',  fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (3, 1, 3)

plt.plot (nn1, g11, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'o', markersize = 5, label = r"$d/a=0.01$")
plt.plot (nn1, g21, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'v', markersize = 5, label = r"$d/a=0.05$")
plt.plot (nn1, g31, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 's', markersize = 5, label = r"$d/a=0.10$")

plt.axvline (q0, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.title  ("$m=1$", fontsize = "15")
plt.xlabel (r'$n\,\epsilon_a$',  fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout()

plt.show()    
#plt.savefig("Fig5.pdf")
