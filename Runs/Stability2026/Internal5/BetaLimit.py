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

infile = open("BetaLimit3.out", "r")

Ec  = []
pc  = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c4 = float(numbers[3])

    Ec.append(c1)
    pc.append(c4)
        
fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim(0., 1.5)
plt.ylim(0., 0.12)

plt.plot (Ea, pa,  color = 'black',  linewidth = 2,   linestyle = 'dotted',  marker = 's', fillstyle = 'none', markersize = 5, label = r"$T_a=0.00$")
plt.plot (Eb, pb,  color = 'black',  linewidth = 2,   linestyle = 'dotted',  marker = '^', fillstyle = 'none', markersize = 5, label = r"$T_a=0.25$")
#plt.plot (Ec, pc,  color = 'black',  linewidth = 2,   linestyle = 'solid')
plt.xlabel (r'$E_a$',               fontsize = "15")
plt.ylabel (r"$p_{0\,{\rm crit}}$", fontsize = "15")
plt.legend (fontsize = 15)

plt.tight_layout ()

plt.show ()    
#plt.savefig("BetaLimit1.pdf")
