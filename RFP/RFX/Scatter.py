import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("beta0_locking_stat_output_Fitzpatrick.txt", "r")

shot    = []
beta0   = []
ip      = []
F       = []
dphi01  = []
minmax0 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    shot.append(c1)
    beta0.append(c2)
    ip.append(c3)
    F.append(c4)
    dphi01.append(c5)
    minmax0.append(c6)

dphi01 = np.multiply(dphi01, 1./math.pi)
    
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot(1, 1, 1)

plt.xlim (0., 0.08)
plt.ylim (-0.25, 0.25)

plt.scatter(beta0, dphi01, 2., marker='s', color = 'black')

#plt.axhline (0., color='black',    linewidth=2., linestyle='dotted')
#plt.axvline (0., color='black',    linewidth=2., linestyle='dotted')

plt.xlabel(r'$\beta_0$', fontsize="20")
plt.ylabel(r'$(\varphi_0-\varphi_1)/\pi$', fontsize="20")
#plt.legend(fontsize="17")

plt.tight_layout()

#plt.show ()
plt.savefig("Figure10_9.pdf")    
