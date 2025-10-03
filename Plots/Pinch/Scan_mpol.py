# Scan_mpol.py

# Plots growth rates versus toroidal mode number

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

infile = open("../../Outputs/Pinch/mscan.out", "r")

nn = []
g1 = []
g2 = []
g3 = []

for line in infile: 

    numbers = line.split()
    q0      = float(numbers[1])
    qa      = float(numbers[3])
    c1      = float(numbers[0])
    c2      = float(numbers[5])
    c3      = float(numbers[6])
    c4      = float(numbers[7])
    nn.append(c1)
    g1.append(c2)
    g2.append(c3)
    g3.append(c4)

fig = plt.figure (figsize = (8., 6.))
fig.canvas.manager.set_window_title (r'Pinch Code: Poloidal Mode Number Scan: n = %3d' % q0)
plt.rc ('xtick', labelsize=15) 
plt.rc ('ytick', labelsize=15)

plt.subplot (1, 1, 1)

plt.plot (nn, g1, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'o', markersize = 5, label = r"$d/a=0.01$")
plt.plot (nn, g2, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'v', markersize = 5, label = r"$d/a=0.05$")
plt.plot (nn, g3, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 's', markersize = 5, label = r"$d/a=0.10$")

plt.axhline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',              fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")
plt.legend (fontsize = "15");

plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()

plt.show()    
#plt.savefig("Figure9_12.pdf")
