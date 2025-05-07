import numpy as np
import matplotlib.pyplot as plt

filename = 'T0infty.out'

data = np.loadtxt(filename, skiprows=1)

x = data[:, 0]
y = data[:, 1]

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.xlim (0., 1.)

plt.plot(x, y, marker='o', linewidth = 2, color = 'black')

#plt.axvline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel('$|\delta|$', fontsize = fontsize)
plt.ylabel('$\delta T_{0\,\infty}$', fontsize = fontsize)

#plt.show()
plt.savefig ("T0infty.pdf")
