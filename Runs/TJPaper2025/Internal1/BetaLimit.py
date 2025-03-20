import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("BetaLimit.out", "r")

bw  = []
pc1 = []
dW1 = []
pc2 = []
dW2 = []

pc  = []
bc  = []

next(infile)
next(infile)
next(infile)

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c2 = float(numbers[1])
    c3 = float(numbers[2])
    c4 = float(numbers[3])
    c5 = float(numbers[4])

    bw.append(c1)
    pc1.append(c2)
    dW1.append(c3)
    pc2.append(c4)
    dW2.append(c5)

b0 = []
for b,p1,w1,p2,w2 in zip(bw,pc1,dW1,pc2,dW2):
    p0 = (p1*w2 - p2*w1) /(w2 - w1)
    pc.append(p2)
    bc.append(0.08*p0)

"""
for b,p in zip(bw, pc):
    print (b, p)
"""
    
fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.plot (bw, bc,  color = 'black',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
plt.axhline (0.0070930, color = 'black', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (1., color = 'black', linewidth = 1.5, linestyle = 'dashed')
plt.xlabel (r'$b_w$',                   fontsize = "15")
plt.ylabel (r"$\beta_{0\,{\rm crit}}$", fontsize = "15")

plt.tight_layout ()

#plt.show ()    
plt.savefig("BetaLimit1.pdf")
