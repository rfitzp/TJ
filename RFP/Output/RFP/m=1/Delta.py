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
    
fig = plt.figure(figsize=(12., 6.))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)
plt.xlim (8., 60.)
plt.ylim (0., 15.)
plt.plot(n0, d0, color='red', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.00$")
plt.plot(n0, d0, 'ko', markersize=2)

plt.plot(n1, d1, color='green', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.01$")
plt.plot(n1, d1, 'ko', markersize=2)

plt.plot(n2, d2, color='blue', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.02$")
plt.plot(n2, d2, 'ko', markersize=2)

plt.plot(n3, d3, color='yellow', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.03$")
plt.plot(n3, d3, 'ko', markersize=2)

plt.plot(n4, d4, color='cyan', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.04$")
plt.plot(n4, d4, 'ko', markersize=2)

plt.plot(n5, d5, color='black', linewidth = 2, linestyle = 'solid', label=r"$\beta_0=0.045$")
plt.plot(n5, d5, 'co', markersize=2)

plt.xlabel(r'$n$', fontsize="20")
plt.ylabel(r'$\Delta^{1,n}$',  fontsize="20")

plt.legend(fontsize="20")

plt.tight_layout()

#plt.show()    

plt.savefig("Figure6.pdf")
