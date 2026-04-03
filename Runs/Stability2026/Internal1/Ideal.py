import math
import numpy as np
import matplotlib.pyplot as plt

infile = open("Ideal.out", "r")

p0  = []
dW  = []
dWp = []
dWv = []

for line in infile: 

    numbers = line.split()
    
    c1 = float(numbers[0])
    c2 = float(numbers[1])
    c3 = float(numbers[2])
    c4 = float(numbers[3])

    p0.append(c1)
    dW.append(c2)
    dWp.append(c3)
    dWv.append(c4)
    
fontsize = 15

fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0., 0.12)
#plt.ylim (0.,  0.08)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.plot (p0, dW,  color = 'black',  linewidth = 2, linestyle = 'solid',   label =  r"$\delta W$")
plt.plot (p0, dWp, color = 'black',  linewidth = 2, linestyle = 'dashed',  label =  r"$\delta W_p$")
plt.plot (p0, dWv, color = 'black',  linewidth = 2, linestyle = 'dotted',  label =  r"$\delta W_v$")

plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')
plt.axvline (0.1104,  color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$p_0$', fontsize = fontsize)
#plt.ylabel (r"$\Delta'$", fontsize = fontsize)
plt.legend (fontsize = 13)

plt.tight_layout ()

plt.show ()   
