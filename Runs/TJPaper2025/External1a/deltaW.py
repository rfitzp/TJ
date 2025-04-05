import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("deltaW.out", "r")

qa = []
qc = []
pc = []
mu = []
ep = []

Wt = []
Wv = []
Wp = [] 

Wt1 = []
Wv1 = []
Wp1 = [] 

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

    c6 = float(numbers[5])
    c7 = float(numbers[6])
    c8 = float(numbers[7])

    Wp.append(c6)
    Wv.append(c7)
    Wt.append(c8)

    c9  = float(numbers[8])
    c10 = float(numbers[9])
    c11 = float(numbers[10])

    Wp1.append(0.8*c9)
    Wv1.append(0.8*c10)
    Wt1.append(0.8*c11)

b0 = []
for p in pc:
    b0.append(0.08*p)

pmin = 0.
pmax = 1.5e-2
ymax = 1.e-6
    
fig = plt.figure (figsize = (8.0, 6.0))
fig.canvas.manager.set_window_title (r'Circular: q_0=1.5: p_sig=1.315, epsilon_a=0.2')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (pmin, pmax)
#plt.ylim (-10., 20000.)

plt.plot (b0, Wp,  color = 'blue',  linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W_p$ (STRIDE)')
plt.plot (b0, Wp1, color = 'blue',  linewidth = 1,  linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5, label = '$\delta W_p$ (TJ)')
plt.plot (b0, Wv,  color = 'red',   linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W_v$ (STRIDE)')
plt.plot (b0, Wv1, color = 'red',   linewidth = 1,  linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5, label = '$\delta W_v$ (TJ)')
plt.plot (b0, Wt,  color = 'black', linewidth = 1,  linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W_t$ (STRIDE)')
plt.plot (b0, Wt1, color = 'black', linewidth = 1,  linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5, label = '$\delta W_t$ (TJ)')
plt.axhline (0.,            color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',    fontsize = "15")
plt.legend (fontsize = "12")

plt.tight_layout ()

#plt.show ()    
plt.savefig("deltaW.pdf")
