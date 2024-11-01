import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Test9.out", "r")

qa = []
qc = []
pc = []
mu = []
ep = []

m1 = []
p1 = []
r1 = []
d1 = []

m2 = []
p2 = []
r2 = []
d2 = []

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
    c9 = float(numbers[8])

    m1.append(c6)
    p1.append(c7)
    r1.append(c8)
    d1.append(c9)

    c10 = int  (numbers[9])
    c11 = float(numbers[10])
    c12 = float(numbers[11])
    c13 = float(numbers[12])

    m2.append(c10)
    p2.append(c11)
    r2.append(c12)
    d2.append(c13)

    c15 = float(numbers[14])
    c16 = float(numbers[15])
    c17 = float(numbers[16])
    c18 = float(numbers[17])

    e11r.append(c15)
    e11i.append(c16)
    E11r.append(c17)
    E11i.append(c18)

    c20 = float(numbers[19])
    c21 = float(numbers[20])
    c22 = float(numbers[21])
    c23 = float(numbers[22])

    e12r.append(c20)
    e12i.append(c21)
    E12r.append(c22)
    E12i.append(c23)

    c25 = float(numbers[24])
    c26 = float(numbers[25])
    c27 = float(numbers[26])
    c28 = float(numbers[27])

    e21r.append(c25)
    e21i.append(-c26)
    E21r.append(c27)
    E21i.append(-c28)

    c30 = float(numbers[29])
    c31 = float(numbers[30])
    c32 = float(numbers[31])
    c33 = float(numbers[32])

    e22r.append(c30)
    e22i.append(c31)
    E22r.append(c32)
    E22i.append(c33)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'q_0=1.1: epsilon_a=0.2: p_0=0.1: H_2=0.0')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (0., 0.2)

plt.plot (mu, e11r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, d1,   color = 'red',   linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TEAR')
plt.plot (mu, E11r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.xlabel (r'$H_2$',        fontsize = "15")
plt.ylabel (r"$Re(E_{11})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (3, 2, 2)

plt.xlim (0., 0.2)

plt.plot (mu, e11i, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, E11i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2$',        fontsize = "15")
plt.ylabel (r"$Im(E_{11})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (3, 2, 3)

plt.xlim (0., 0.2)

plt.plot (mu, e12r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, e21r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)
plt.plot (mu, E12r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')
plt.plot (mu, E21r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)

plt.xlabel (r'$H_2$',                    fontsize = "15")
plt.ylabel (r"$Re(E_{12}), Re(E_{21})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (3, 2, 4)

plt.xlim (0., 0.2)

plt.plot (mu, e12i, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, e21i, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)
plt.plot (mu, E12i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')
plt.plot (mu, E21i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2$',                    fontsize = "15")
plt.ylabel (r"$Im(E_{12}), -Im(E_{21})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (3, 2, 5)

plt.xlim (0., 0.2)

plt.plot (mu, e22r, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, d2,   color = 'red',   linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TEAR')
plt.plot (mu, E22r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.xlabel (r'$H_2$',        fontsize = "15")
plt.ylabel (r"$Re(E_{22})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (3, 2, 6)

plt.xlim (0., 0.2)

plt.plot (mu, e22i, color = 'blue',  linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'TJ')
plt.plot (mu, E22i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = 'STRIDE')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2$',        fontsize = "15")
plt.ylabel (r"$Im(E_{22})$", fontsize = "15")
plt.legend (fontsize = "12")

plt.tight_layout ()

plt.show ()    
#plt.savefig("Test3.pdf")
