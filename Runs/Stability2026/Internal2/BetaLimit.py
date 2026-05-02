import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("BetaLimit.out", "r")

q = []
r = []
s = []
p = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c2 = float(numbers[1])
    c3 = float(numbers[2])
    c4 = float(numbers[3])

    q.append(c1)
    r.append(c2)
    s.append(c3)
    p.append(c4)
    
fontsize = 15

fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (1.0, 1.4)
plt.ylim (0.105, 0.125)

plt.plot (q, p,  color = 'black',  linewidth = 2, linestyle = 'solid')

#plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axhline (1.1046e-01 ,  color = 'black', linewidth = 1, linestyle = 'dashed')

plt.xlabel (r'$b_w$', fontsize = fontsize)
plt.ylabel (r"$p_{0\,cirit}$", fontsize = fontsize)

plt.tight_layout ()

plt.show ()   
