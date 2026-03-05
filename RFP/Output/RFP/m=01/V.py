import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("beta000.out", "r")

n0 = []
k0 = []
x0 = []
y0 = []
z0 = []
j0 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n0.append(c1)
    k0.append(c2)
    x0.append(c3)
    y0.append(c4)
    z0.append(c5)
    j0.append(c6)

infile = open("beta001.out", "r")

n1 = []
k1 = []
x1 = []
y1 = []
z1 = []
j1 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n1.append(c1)
    k1.append(c2)
    x1.append(c3)
    y1.append(c4)
    z1.append(c5)
    j1.append(c6)

infile = open("beta002.out", "r")

n2 = []
k2 = []
x2 = []
y2 = []
z2 = []
j2 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n2.append(c1)
    k2.append(c2)
    x2.append(c3)
    y2.append(c4)
    z2.append(c5)
    j2.append(c6)        

infile = open("beta003.out", "r")

n3 = []
k3 = []
x3 = []
y3 = []
z3 = []
j3 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n3.append(c1)
    k3.append(c2)
    x3.append(c3)
    y3.append(c4)
    z3.append(c5)
    j3.append(c6)        

infile = open("beta004.out", "r")

n4 = []
k4 = []
x4 = []
y4 = []
z4 = []
j4 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n4.append(c1)
    k4.append(c2)
    x4.append(c3)
    y4.append(c4)
    z4.append(c5)
    j4.append(c6)

infile = open("beta0045.out", "r")

n5 = []
k5 = []
x5 = []
y5 = []
z5 = []
j5 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    n5.append(c1)
    k5.append(c2)
    x5.append(c3)
    y5.append(c4)
    z5.append(c5)
    j5.append(c6)      
      
def V0(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k0, x0, y0, z0, j0):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.75

def V1(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k1, x1, y1, z1, j1):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.67

def V2(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k2, x2, y2, z2, j2):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.64

def V3(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k3, x3, y3, z3, j3):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.53

def V4(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k4, x4, y4, z4, j4):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.6

def V5(phi):
    sum = 0.
    norm = 0.
    for k, x, y, z, j in zip(k5, x5, y5, z5, j5):
        sum  -= x*x*y*y*z*z*j * np.sin(k*phi*math.pi)
        norm += x*x*y*y*z*z*j

    return sum/abs(norm)/0.95

fig = plt.figure(figsize=(12., 6.))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

phi = np.linspace(-1., 1., 1000)

plt.subplot(1, 1, 1)
plt.xlim (-1., 1.)
plt.ylim (-1.0, 1.0)
plt.plot(phi, V0(phi), color = 'red', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.00$")
plt.plot(phi, V1(phi), color = 'green', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.03$")
plt.plot(phi, V2(phi), color = 'blue', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.06$")
plt.plot(phi, V3(phi), color = 'yellow', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.09$")
plt.plot(phi, V4(phi), color = 'cyan', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.04$")
plt.plot(phi, V5(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.045$")
plt.axhline(y=0., color = 'black', linestyle = 'dotted', linewidth = '1.')
plt.axvline(x=0., color = 'black', linestyle = 'dotted', linewidth = '1.')

plt.xlabel(r'$(\phi_1-\phi_0)/\pi$', fontsize="20")
plt.ylabel(r'$V(a.u.)$',  fontsize="20")

plt.legend(fontsize="15")

plt.tight_layout()

plt.show()

#plt.savefig("Figure8.pdf")
