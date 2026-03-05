import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("brm0_Delta1_0.100000_19596.0.txt", "r")

phix    = []
brx     = []
Deltax  = []
brhx    = []
Deltahx = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    phix.append(c1)
    brx.append(c2)
    Deltax.append(c3)
    brhx.append(c4)
    Deltahx.append(c5)

px = np.multiply(phix, 1./math.pi)    

infile = open ("brm0_Delta1_0.100000_29509.0.txt", "r")

phi    = []
br     = []
Delta  = []
brh    = []
Deltah = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    phi.append(c1)
    br.append(c2)
    Delta.append(c3)
    brh.append(c4)
    Deltah.append(c5)

p = np.multiply(phi, 1./math.pi)
p = np.subtract(p, 1.)

brh1     = np.roll(brh, 180)
Deltah1  = np.roll(Deltah, 180)

fig = plt.figure (figsize=(12.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot(1, 2, 1)

plt.xlim (0., 2.)
plt.ylim (-0.65, 1.05)

plt.plot (px, brhx,    color='black',  linewidth = 2, linestyle = 'solid',  label = r'$b_r^{\,m=0}$')
plt.plot (px, Deltahx, color='black',  linewidth = 2, linestyle = 'dashed', label = r'$\Delta^{m=1}$')

plt.axhline (0., color='black',    linewidth=2., linestyle='dotted')
#plt.axvline (0., color='black',    linewidth=2., linestyle='dotted')

#plt.xlabel(r'$\hat{x}$', fontsize="20")
plt.title(r'19596', fontsize='20')
plt.xlabel(r'$\varphi/\pi$', fontsize="20")
plt.legend(fontsize="17")

plt.subplot(1, 2, 2)

plt.xlim (-1., 1.)
plt.ylim (-0.65, 1.05)

plt.plot (p, brh1,    color='black', linewidth = 2, linestyle = 'solid',  label = r'$b_r^{\,m=0}$')
plt.plot (p, Deltah1, color='black', linewidth = 2, linestyle = 'dashed', label = r'$\Delta^{m=1}$')

plt.axhline (0., color='black',    linewidth=2., linestyle='dotted')
#plt.axvline (0., color='black',    linewidth=2., linestyle='dotted')

#plt.xlabel(r'$\hat{x}$', fontsize="20")
plt.title(r'29509', fontsize='20')
plt.xlabel(r'$\varphi/\pi$', fontsize="20")
plt.legend(fontsize="17")

plt.tight_layout(pad=0.5)

#plt.show ()
plt.savefig ("Figure10_8.pdf")
