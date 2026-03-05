import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("beta01.out", "r")

n0 = []
d0 = []
p0 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    n0.append(c1)
    d0.append(c2)
    p0.append(c3)

infile = open("beta03.out", "r")

n1 = []
d1 = []
p1 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    n1.append(c1)
    d1.append(c2)
    p1.append(c3)

infile = open("beta06.out", "r")

n2 = []
d2 = []
p2 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    n2.append(c1)
    d2.append(c2)
    p2.append(c3)

infile = open("beta09.out", "r")

n3 = []
d3 = []
p3 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    n3.append(c1)
    d3.append(c2)
    p3.append(c3)

infile = open("beta12.out", "r")

n4 = []
d4 = []
p4 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    n4.append(c1)
    d4.append(c2)
    p4.append(c3)

infile = open("beta15.out", "r")

n5 = []
d5 = []
p5 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    n5.append(c1)
    d5.append(c2)
    p5.append(c3)

def br0(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n0, d0, p0):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm

def br1(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n1, d1, p1):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

def br2(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n2, d2, p2):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

"""
def br3(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n3, d3, p3):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

def br4(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n4, d4, p4):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm

def br5(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n5, d5, p5):
        sum += d*d * p * np.cos(n*phi*math.pi)
        norm += d*d * p

    return sum/norm
"""

def cr0(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n0, d0, p0):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm

def cr1(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n1, d1, p1):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

def cr2(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n2, d2, d3):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

"""
def cr3(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n3, d3, p3):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  

def cr4(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n4, d4, p4):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm

def cr5(phi):
    sum = 0.
    norm = 0.
    for n, d, p in zip(n5, d5, p5):
        sum -= d*d * p * np.sin(n*phi*math.pi)
        norm += d*d * p

    return sum/norm  
"""
fig = plt.figure(figsize=(8., 10.))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

phi = np.linspace(-1., 1., 1000)

plt.subplot(2, 1, 1)
plt.xlim (-1., 1.)
plt.ylim (-1.1, 1.1)
plt.plot(phi, br0(phi), color = 'black', linewidth = 2, linestyle = 'dotted', label=r"$\beta_0=0.01$")
plt.plot(phi, br1(phi), color = 'black', linewidth = 2, linestyle = 'dashed', label=r"$\beta_0=0.03$")
plt.plot(phi, br2(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.06$")
#plt.plot(phi, br3(phi), color = 'yellow', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.09$")
#plt.plot(phi, br4(phi), color = 'cyan', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.12$")
#plt.plot(phi, br5(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.15$")
plt.axhline(y=0., color = 'black', linestyle = 'dotted', linewidth = '1.')
plt.axvline(x=0., color = 'black', linestyle = 'dotted', linewidth = '1.')

plt.xlabel(r'$(\varphi-\varphi_1)/\pi$', fontsize="20")
plt.ylabel(r'$C^{\,m=1}(a.u.)$',  fontsize="20")

plt.legend(fontsize="15")

plt.subplot(2, 1, 2)
plt.xlim (-1., 1.)
plt.ylim (-1.1, 1.1)
plt.plot(phi, cr0(phi), color = 'black', linewidth = 2, linestyle = 'dotted', label=r"$\beta_0=0.01$")
plt.plot(phi, cr1(phi), color = 'black', linewidth = 2, linestyle = 'dashed', label=r"$\beta_0=0.03$")
plt.plot(phi, cr2(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.06$")
#plt.plot(phi, cr3(phi), color = 'yellow', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.09$")
#plt.plot(phi, cr4(phi), color = 'cyan', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.12$")
#plt.plot(phi, cr5(phi), color = 'black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.15$")
plt.axhline(y=0., color = 'black', linestyle = 'dotted', linewidth = '1.')
plt.axvline(x=0., color = 'black', linestyle = 'dotted', linewidth = '1.')

plt.xlabel(r'$(\varphi-\varphi_1)/\pi$', fontsize="20")
plt.ylabel(r'$S^{\,m=1}(a.u.)$',  fontsize="20")

plt.legend(fontsize="15")

plt.tight_layout()

#plt.show()

plt.savefig("Figure10_7.pdf")
