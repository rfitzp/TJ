import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Test3.out", "r")

qa = []
pn = []
rs = []
d1 = []
d2 = []
d3 = []

next(infile)
next(infile)

for line in infile: 

    numbers = line.split()
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    qa.append(c1)
    pn.append(c2)
    rs.append(c3)
    d1.append(c4)
    d2.append(c5)
    d3.append(c6)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Circular: Zero pressure: q_0=1.1: epsilon_a=0.2')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (2.2, 3.)

plt.plot (qa, d1, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = 'TJ')
plt.plot (qa, d2, color = 'red',   linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = 'TEAR')
plt.plot (qa, d3, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = 'STRIDE')

plt.xlabel (r'$q_a$',          fontsize = "15")
plt.ylabel (r"$r_s\,\Delta'$", fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
#plt.savefig("Test3.pdf")
