import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Test4.out", "r")

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
E11r = []
E11i = []

e12r = []
e12i = []
E12r = []
E12i = []

e21r = []
e21i = []
E21r = []
E21i = []

e22r = []
e22i = []
E22r = []
E22i = []

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
    c17 = float(numbers[14])
    c18 = float(numbers[15])

    e11r.append(c15)
    e11i.append(c16)
    E11r.append(c17)
    E11i.append(c18)

    c20 = float(numbers[17])
    c21 = float(numbers[18])
    c22 = float(numbers[19])
    c23 = float(numbers[20])

    e12r.append(c20)
    e12i.append(c21)
    E12r.append(c22)
    E12i.append(c23)

    c25 = float(numbers[22])
    c26 = float(numbers[23])
    c27 = float(numbers[24])
    c28 = float(numbers[25])

    e21r.append(c25)
    e21i.append(-c26)
    E21r.append(c27)
    E21i.append(-c28)

    c30 = float(numbers[27])
    c31 = float(numbers[28])
    c32 = float(numbers[29])
    c33 = float(numbers[30])

    e22r.append(c30)
    e22i.append(c31)
    E22r.append(c32)
    E22i.append(c33)

b0 = []
b1 = []
for p in pc:
    b0.append(0.08*p)
    b1.append(0.08*p*1.08)

pmin = 0.
pmax = 2.5e-2
ymax = 1.e-6
    
fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'Circular: q_0=0.8: epsilon_a=0.2')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 1, 1)

plt.xlim (pmin, pmax)
plt.ylim (0., 100.)

plt.plot (b1, e11r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (b0, E11r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.xlabel (r'$\beta_0$', fontsize = "17")
plt.ylabel (r"$E_{11}$",  fontsize = "17")
plt.legend (fontsize = "12")

plt.subplot (3, 1, 2)

plt.xlim (pmin, pmax)
plt.ylim (-80., 0.)

plt.plot (b1, e12r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (b1, e21r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)
plt.plot (b0, E12r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')
plt.plot (b0, E21r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)

plt.xlabel (r'$\beta_0$',        fontsize = "17")
plt.ylabel (r"$E_{11}, E_{21}$", fontsize = "17")
plt.legend (fontsize = "12", loc = "lower left")

plt.subplot (3, 1, 3)

plt.xlim (pmin, pmax)
plt.ylim (-5., 50.)

plt.plot (b1, e22r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (b0, E22r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.xlabel (r'$\beta_0$', fontsize = "17")
plt.ylabel (r"$E_{22}$",  fontsize = "17")
plt.legend (fontsize = "12")

plt.tight_layout ()

#plt.show ()    
plt.savefig("Test4.pdf")
