import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Overlap.out", "r")

rr1  = []
t1   = []
tp1  = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    rr1.append(c1)
    t1.append(c2)
    tp1.append(c3)
     
infile = open("Rational.out", "r")

for line in infile: 

    numbers = line.split() 
    rs1     = float(numbers[0])
    rs2     = float(numbers[1])
    rs3     = float(numbers[2])

infile = open("Integral.out", "r")

rr  = []
t   = []
tp  = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    rr.append(c1)
    t.append(c2)
    tp.append(c3)
     
fig = plt.figure(figsize=(12.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(2, 1, 1)

plt.xlim(0., 1.)

plt.plot(rr1, tp1, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr1, t1, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs1, color='red', linewidth=2., linestyle='dotted')
plt.axvline(rs2, color='red', linewidth=2., linestyle='dotted')
plt.axvline(rs3, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r'$\hat{\tau}$',  fontsize="20")

plt.subplot(2, 1, 2)

plt.xlim(0., 1.)

plt.plot(rr, tp, color='blue', linewidth = 2, linestyle = 'solid')
plt.plot(rr, t, color='blue', linewidth = 2, linestyle = 'dotted')
plt.axvline(rs1, color='red', linewidth=2., linestyle='dotted')
plt.axvline(rs2, color='red', linewidth=2., linestyle='dotted')
plt.axvline(rs3, color='red', linewidth=2., linestyle='dotted')
plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

plt.xlabel(r'$\hat{r}$', fontsize="20")
plt.ylabel(r"$\int_0^{\hat{r}}\hat{\tau}(\hat{r}')\,d\hat{r}'$",  fontsize="20")

plt.tight_layout()

plt.show()    

#plt.savefig("Figure5.pdf")
