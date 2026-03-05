import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Control.out", "r")

cc = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    cc.append(c1)

infile = open("Eigenfunction.out", "r")

rr  = []
p1  = []
pp1 = []
p2  = []
pp2 = []
p3  = []
pp3 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    c7      = float(numbers[6])
    rr.append(c1)
    p1.append(c2)
    pp1.append(c3)
    p2.append(c4)
    pp2.append(c5)
    p3.append(c6)
    pp3.append(c7)

infile = open("EigenStepped.out", "r")

xrr  = []
xp1  = []
xpp1 = []
xp2  = []
xpp2 = []
xp3  = []
xpp3 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    c7      = float(numbers[6])
    xrr.append(c1)
    xp1.append(c2)
    xpp1.append(c3)
    xp2.append(c4)
    xpp2.append(c5)
    xp3.append(c6)
    xpp3.append(c7)
    
infile = open("Rational.out", "r")

for line in infile: 

    numbers = line.split() 
    rs1     = float(numbers[0])
    rs2     = float(numbers[1])
    rs3     = float(numbers[2])

infile = open("Mode.out", "r")

for line in infile: 

    numbers = line.split() 
    m1     = int(numbers[0])
    n1     = int(numbers[1])
    m2     = int(numbers[2])
    n2     = int(numbers[3])
    m3     = int(numbers[4])
    n3     = int(numbers[5])
                      
fig = plt.figure(figsize=(12.0, 10.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(3, 2, 1)

plt.xlim(0., 1.)

plt.plot(xrr, xp1, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, p1, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs1, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')
if m1 == 1 and n1 == 1:
    pass
else:
    plt.axhline(1., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r'$\hat{\psi}^{\,m_1,n_1}$',  fontsize="20")

plt.subplot(3, 2, 2)

plt.xlim(0., 1.)

plt.plot(xrr, xpp1, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, pp1, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs1, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r"$(\hat{\psi}^{\,m_1,n_1})'$",  fontsize="20")

plt.subplot(3, 2, 3)

plt.xlim(0., 1.)

plt.plot(xrr, xp2, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, p2, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs2, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')
if m2 == 1 and n2 == 1:
    pass
else:
    plt.axhline(1., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r'$\hat{\psi}^{\,m_2,n_2}$',  fontsize="20")

plt.subplot(3, 2, 4)

plt.xlim(0., 1.)

plt.plot(xrr, xpp2, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, pp2, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs2, color='black', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r"$(\hat{\psi}^{\,m_2,n_2})'$",  fontsize="20")

plt.subplot(3, 2, 5)

plt.xlim(0., 1.)

plt.plot(xrr, xp3, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, p3, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs3, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')
if m3 == 1 and n3 == 1:
    pass
else:
    plt.axhline(1., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r'$\hat{\psi}^{\,m_3,n_3}$',  fontsize="20")

plt.subplot(3, 2, 6)

plt.xlim(0., 1.)

plt.plot(xrr, xpp3, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, pp3, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs3, color='red', linewidth=2.0, linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r"$(\hat{\psi}^{\,m_3,n_3})'$",  fontsize="20")

plt.tight_layout()

#plt.show()    

plt.savefig("Figure2.pdf")
