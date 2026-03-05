import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Control.out", "r")

cc = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    cc.append(c1)

infile = open("Stepped.out", "r")

rr = []
ss = []
pp = []
bp = []
bt = []
qq = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    rr.append(c1)
    ss.append(c2)
    pp.append(c3)
    bp.append(c4)
    bt.append(c5)
    qq.append(c6)

infile = open("Equilibrium.out", "r")

rrx = []
ssx = []
ppx = []
bpx = []
btx = []
qqx = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    rrx.append(c1)
    ssx.append(c2)
    ppx.append(c3)
    bpx.append(c4)
    btx.append(c5)
    qqx.append(c6)

infile = open("Zero.out", "r")

rry = []
ssy = []
ppy = []
bpy = []
bty = []
qqy = []
pq  = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    c7      = float(numbers[6])
    rry.append(c1)
    ssy.append(c2)
    ppy.append(c3)
    bpy.append(c4)
    bty.append(c5)
    qqy.append(c6)
    pq.append(c7) 

infile = open("Rational.out", "r")

for line in infile: 

    numbers = line.split() 
    rs1      = float(numbers[0])
    rs2      = float(numbers[1])
    rs3      = float(numbers[2])
        
fig = plt.figure(figsize=(12., 8.))
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 

plt.subplot(2, 2, 1)

plt.xlim(0., 1.)
plt.ylim(0., 1.05*ss[0])

plt.plot(rr, ss, color='black', linewidth = 2, linestyle = 'solid')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\bar{r}$', fontsize="20")
plt.ylabel(r'$\sigma$',  fontsize="20")

plt.subplot(2, 2, 2)

plt.xlim(0., 1.)
plt.ylim(0., 1.05*ppy[0])

#plt.plot(rrx, ppx, color='black', linewidth = 2, linestyle = 'dotted')
plt.plot(rr, ppy, color='black', linewidth = 2, linestyle = 'solid')
#plt.plot(rr, pq, color='green', linewidth = 2, linestyle = 'solid')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\bar{r}$',  fontsize="20")
plt.ylabel(r'$\bar{P}$',  fontsize="20")

plt.subplot(2, 2, 3)

plt.xlim(0., 1.)

if bp[-1] > 0.:
    plt.ylim(0., 1.05)
else:
    plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')

plt.plot(rrx, bpx, color='black', linewidth = 2, linestyle = 'dotted')
plt.plot(rrx, btx, color='black', linewidth = 2, linestyle = 'dotted')
plt.plot(rr, bp, color='black', linewidth = 2, linestyle = 'solid', label = r'$\bar{B}_\varphi$')
plt.plot(rr, bt, color='black', linewidth = 2, linestyle = 'dashed', label = r'$\bar{B}_\theta$')
plt.plot(rry, bpy, color='black', linewidth = 2, linestyle = 'dotted')
plt.plot(rry, bty, color='black', linewidth = 2, linestyle = 'dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\bar{r}$', fontsize="20")
plt.legend(fontsize="20")

plt.subplot(2, 2, 4)

plt.xlim(0., 1.)

if bp[-1] > 0.:
    pass
else:
    plt.axhline(0., color='black', linewidth=1.5, linestyle='dotted')

plt.plot(rr, qq, color='black', linewidth = 2, linestyle = 'solid')
#plt.plot(rrx, qqx, color='black', linewidth = 2, linestyle = 'dotted')
#plt.axvline(rs1, color='red', linewidth=2., linestyle='dotted')
#plt.axvline(rs2, color='red', linewidth=2., linestyle='dotted')
#plt.axvline(rs3, color='red', linewidth=2., linestyle='dotted')

if (len(cc) < 50):
    for c in cc:
        plt.axvline(c, color='black', linewidth=1., linestyle='dotted')

plt.xlabel(r'$\bar{r}$', fontsize="20")
plt.ylabel(r'$q$',  fontsize="20")

plt.tight_layout()

#plt.show()    
plt.savefig("Figure1.pdf")
