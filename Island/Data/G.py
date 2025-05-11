import numpy as np
import matplotlib.pyplot as plt

filename = 'delta.out'

data = np.loadtxt(filename)

delta = data[:, 0]
gruth = data[:, 1]
gboot = data[:, 3]/6.

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.xlim (0., 1.)
plt.ylim (0., 1.2)

plt.plot (delta, gruth,     color = 'red',     linewidth = 2, linestyle = 'solid', label = r"$G_{ruth}$")
plt.plot (delta, gboot,     color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$G_{boot}/6$")

#plt.axvline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel('$|\delta|$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

#plt.show()
plt.savefig ("G.pdf")
