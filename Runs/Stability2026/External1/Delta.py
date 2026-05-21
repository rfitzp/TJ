import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

infile = open("TJ.out", "r")

p2 = []
D1 = []
D2 = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c2 = float(numbers[5])
    c3 = float(numbers[7])
    c4 = float(numbers[8])
    c5 = float(numbers[9])

    p2.append(c1)
    D1.append(math.log10(c2))
    D2.append(math.log10(c3))
 
fontsize = 15

fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0., 0.18)
#plt.ylim (0.,  0.08)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.plot (p2, D1, color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"$E_{11}$")
plt.plot (p2, D2, color = 'black',  linewidth = 2, linestyle = 'dashed',  label =  r"$E_{22}$")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axvline (0.1674,  color = 'black', linewidth = 1, linestyle = 'dashed')

plt.xlabel (r'$p_0$', fontsize = fontsize)
plt.ylabel (r"$\log_{10}(E_{kk})$", fontsize = fontsize)
plt.legend (fontsize = 15)

plt.tight_layout ()

plt.show ()    
#plt.savefig ("Figure13_8.pdf")
