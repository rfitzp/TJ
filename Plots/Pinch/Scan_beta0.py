# Scan_beta0.py

# Plots growth rates versus central beta

import math
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    gmax = float (sys.argv[1])
if len(sys.argv) > 2:
    gmin = float (sys.argv[1])
    gmax = float (sys.argv[2])

infile = open("../../Outputs/Pinch/b0scan.out", "r")

bb = []
g1 = []
g2 = []
g3 = []

for line in infile: 

    numbers = line.split()
    mm      = float(numbers[0])
    nn      = float(numbers[1])
    c1      = float(numbers[5])
    c2      = float(numbers[6])
    c3      = float(numbers[7])
    c4      = float(numbers[8])
    bb.append(c1)
    g1.append(c2)
    g2.append(c3)
    g3.append(c4)

fig = plt.figure (figsize = (8., 6.))
fig.canvas.manager.set_window_title (r'Pinch Code: Central beta Scan: m = %2d, n = %2d' % (mm, nn))
plt.rc ('xtick', labelsize=15) 
plt.rc ('ytick', labelsize=15)

plt.subplot (1, 1, 1)

if (len(sys.argv) > 1):
    plt.ylim (0., gmax)
if (len(sys.argv) > 2):
    plt.ylim (gmin, gmax)

plt.plot (bb, g1, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'o', markersize = 5, label = r"$d/a=0.01$")
plt.plot (bb, g2, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 'v', markersize = 5, label = r"$d/a=0.05$")
plt.plot (bb, g3, color = 'black', linewidth = 0.5, linestyle = 'dotted', markerfacecolor = 'none', marker = 's', markersize = 5, label = r"$d/a=0.10$")

plt.axhline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout()

plt.show()    
#plt.savefig("Fig5.pdf")
