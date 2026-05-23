import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("BetaLimit1.out", "r")

Ea  = []
pa  = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c4 = float(numbers[3])

    Ea.append(c1)
    pa.append(c4)

infile = open("BetaLimit2.out", "r")

Eb  = []
pb  = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c4 = float(numbers[3])

    Eb.append(c1)
    pb.append(c4)

infile = open("BetaLimit.out", "r")

Ec  = []
pc  = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c4 = float(numbers[3])

    Ec.append(c1)
    pc.append(c4)
    
fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 2, 1)

plt.xlim(0., 1.0)
#plt.ylim(0., 0.10)

plt.plot (Ea, pa,  color = 'black',  linewidth = 2,   linestyle = 'dotted',  marker = 's', fillstyle = 'none', markersize = 5)
plt.xlabel (r'$E_a$',               fontsize = "20")
plt.ylabel (r"$p_{0\,{\rm crit}}$", fontsize = "20")

plt.subplot (1, 2, 2)

plt.xlim(-0.25, 0.3)

plt.plot (Eb, pb,  color = 'black',  linewidth = 2,   linestyle = 'dotted',  marker = 's', fillstyle = 'none', markersize = 5)
plt.plot (Ec, pc,  color = 'black',  linewidth = 2,   linestyle = 'dotted',  marker = 's', fillstyle = 'none', markersize = 5)
plt.axvline (0., linewidth = 2,   linestyle = 'dotted', color = 'black')
plt.xlabel (r'$T_a$',               fontsize = "20")
plt.ylabel (r"$p_{0\,{\rm crit}}$", fontsize = "20")

plt.tight_layout ()

#plt.show ()    
plt.savefig("Figure13_18.pdf")
