import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("deltaW.out", "r")

qa = []
qc = []
pc = []
mu = []
ep = []

W  = []
Wv = []
Wp = [] 

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

    W.append(c6)
    Wp.append(c7)
    Wv.append(c8)

b0 = []
for p in pc:
    b0.append(0.08*p)

pmin = 0.
pmax = 7.5e-3
ymax = 1.e-6
    
fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (pmin, pmax)
#plt.ylim (-10., 20000.)

plt.plot    (b0, W,  color = 'black',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W$')
plt.plot    (b0, Wp, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W_p$')
plt.plot    (b0, Wv, color = 'red',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$\delta W_v$')
plt.axvline (0.088519*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,            color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$\delta W, \delta W_p, \delta W_v$", fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
#plt.savefig("Test4.pdf")
