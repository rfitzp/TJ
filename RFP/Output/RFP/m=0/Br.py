import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("beta000.out", "r")

n0 = []
d0 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n0.append(c1)
    d0.append(c2)

infile = open("beta001.out", "r")

n1 = []
d1 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n1.append(c1)
    d1.append(c2)    

infile = open("beta002.out", "r")

n2 = []
d2 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n2.append(c1)
    d2.append(c2)

infile = open("beta003.out", "r")

n3 = []
d3 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n3.append(c1)
    d3.append(c2)

infile = open("beta004.out", "r")

n4 = []
d4 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n4.append(c1)
    d4.append(c2)

infile = open("beta0045.out", "r")

n5 = []
d5 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n5.append(c1)
    d5.append(c2) 

def br0(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n0, d0):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm

def br1(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n1, d1):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm  

def br2(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n2, d2):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm  

def br3(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n3, d3):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm  

def br4(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n4, d4):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm

def br5(phi):
    sum = 0.
    norm = 0.
    for n, d in zip(n5, d5):
        sum += d*d * np.cos(n*phi*math.pi)
        norm += d*d

    return sum/norm     

fig = plt.figure(figsize=(12., 6.))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

phi = np.linspace(-1., 1., 1000)

plt.subplot(1, 1, 1)
plt.xlim (-1., 1.)
plt.ylim (-0.5, 1.05)
plt.plot(phi, br0(phi), color = 'red', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.00$")
plt.plot(phi, br1(phi), color = 'green', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.01$")
plt.plot(phi, br2(phi), color = 'blue', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.02$")
plt.plot(phi, br3(phi), color = 'yellow', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.03$")
plt.plot(phi, br4(phi), color = 'cyan', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.04$")
plt.plot(phi, br5(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.045$")
plt.axhline(y=0., color = 'black', linestyle = 'dotted', linewidth = '1.')
plt.axvline(x=0., color = 'black', linestyle = 'dotted', linewidth = '1.')

plt.xlabel(r'$(\phi-\phi_0)/\pi$', fontsize="20")
plt.ylabel(r'$\hat{b}_r^{m=0}(a.u.)$',  fontsize="20")

plt.legend(fontsize="20")

plt.tight_layout()

#plt.show()

plt.savefig("Figure4.pdf")
