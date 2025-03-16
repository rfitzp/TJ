import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Ematrix.out", "r")

qa = []
qc = []
pc = []
mu = []
ep = []

m1 = []
p1 = []
r1 = []

m2 = []
p2 = []
r2 = []

e11r = []
e11i = []

e12r = []
e12i = []

e21r = []
e21i = []

e22r = []
e22i = []

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

    c6 = int  (numbers[5])
    c7 = float(numbers[6])
    c8 = float(numbers[7])

    m1.append(c6)
    p1.append(c7)
    r1.append(c8)

    c10 = int  (numbers[8])
    c11 = float(numbers[9])
    c12 = float(numbers[10])

    m2.append(c10)
    p2.append(c11)
    r2.append(c12)

    c15 = float(numbers[12])
    c16 = float(numbers[13])

    e11r.append(c15)
    e11i.append(c16)

    c20 = float(numbers[15])
    c21 = float(numbers[16])

    e12r.append(c20)
    e12i.append(c21)

    c25 = float(numbers[18])
    c26 = float(numbers[19])

    e21r.append(c25)
    e21i.append(-c26)

    c30 = float(numbers[21])
    c31 = float(numbers[22])

    e22r.append(c30)
    e22i.append(c31)

b0 = []
for p in pc:
    b0.append(0.08*p)

pmin = 0.
pmax = 7.5e-3
ymax = 1.e-6
    
fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 1, 1)

plt.xlim (pmin, pmax)
plt.ylim (-10., 20000.)

plt.plot    (b0, e11r, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
plt.axvline (0.0886*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$Re(E_{11})$", fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim (pmin, pmax)
plt.ylim (10., -3000.)

plt.plot (b0, e12r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = '$E_{12}$')
plt.plot (b0, e21r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5, label = '$E_{21}$')
plt.axvline (0.088519*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',                fontsize = "15")
plt.ylabel (r"$Re(E_{12}), Re(E_{21})$", fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (3, 1, 3)

plt.xlim (pmin, pmax)
plt.ylim (-10., 500.)

plt.plot (b0, e22r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
plt.axvline (0.088519*0.08,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\beta_0$',        fontsize = "15")
plt.ylabel (r"$Re(E_{22})$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
#plt.savefig("Test4.pdf")
