import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

infile = open("Layer.out", "r")

p2 = []
g1 = []
w1 = []
f1 = []
g2 = []
w2 = []
f2 = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c2 = float(numbers[1])
    c3 = float(numbers[2])
    c4 = float(numbers[3])
    c5 = float(numbers[4])
    c6 = float(numbers[5])
    c7 = float(numbers[6])

    p2.append(c1)
    g1.append(c2)
    w1.append(c3)
    f1.append(c4)
    g2.append(c5)
    w2.append(c6)
    f2.append(c7)
 
fontsize = 15

fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (3, 1, 1)

plt.xlim (0., 0.12)

plt.plot (p2, g1, color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"m=1")
plt.plot (p2, g2, color = 'black',  linewidth = 2, linestyle = 'dashed',  label =  r"m=2")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axvline (0.1104,  color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$p_0$',        fontsize = fontsize)
plt.ylabel (r"$\gamma(kHz)$", fontsize = fontsize)
plt.legend (fontsize = 13)

plt.subplot (3, 1, 2)

plt.xlim (0., 0.12)
plt.plot (p2, w1, color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"m=1")
plt.plot (p2, w2, color = 'black',  linewidth = 2, linestyle = 'dashed',  label =  r"m=2")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axvline (0.1104,  color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$p_0$',        fontsize = fontsize)
plt.ylabel (r"$\omega(kHz)$", fontsize = fontsize)
plt.legend (fontsize = 13)

plt.subplot (3, 1, 3)

plt.xlim (0., 0.12)
plt.plot (p2, f1, color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"m=1")
plt.plot (p2, f2, color = 'black',  linewidth = 2, linestyle = 'dashed',  label =  r"m=2")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axvline (0.1104,  color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$p_0$',  fontsize = fontsize)
plt.ylabel (r"$f$",    fontsize = fontsize)
plt.legend (fontsize = 13)

plt.tight_layout ()

plt.show ()    
#plt.savefig ("Figure9_16.pdf")
