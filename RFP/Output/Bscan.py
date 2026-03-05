import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Bscan.out", "r")

rr = []
ss = []
pp = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    rr.append(c1)
    ss.append(c2)
    pp.append(c3)
        
fig = plt.figure(figsize=(12., 6.))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)
#plt.xlim (0., rr[-1])
#plt.ylim (0., 1.05*pp[-1])
plt.plot(rr, ss, color='blue', linewidth = 2, linestyle = 'dotted')
plt.plot(rr, pp, color='blue', linewidth = 2, linestyle = 'solid')

plt.xlabel(r'$\beta_0$', fontsize="20")
plt.ylabel(r'$\int_0^1 \hat{\tau}\,d\hat{r}$',  fontsize="20")

plt.tight_layout()

plt.show()    
