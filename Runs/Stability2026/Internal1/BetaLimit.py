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

plt.xlim (0.9, 1.0)
#plt.ylim (0.,  0.08)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.plot (q, p,  color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"$p_{0\,crit}$")
plt.plot (q, r,  color = 'black',  linewidth = 2, linestyle = 'dashed',   label =  r"$r_s$")
plt.plot (q, s,  color = 'black',  linewidth = 2, linestyle = 'dotted',   label =  r"$s_s$")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
#plt.axvline (0.1104,  color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$q_0$', fontsize = fontsize)
#plt.ylabel (r"$\Delta'$", fontsize = fontsize)
plt.legend (fontsize = 13)

plt.tight_layout ()

plt.show ()   
