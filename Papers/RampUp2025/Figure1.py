import matplotlib.pyplot as plt
import numpy as np

font = 17
fig  = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font)

theta = np.linspace(0, 2*np.pi, 400)

a = 0.32
b = 0.32*1.8
x = 1. + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 4, color = "black")

a = 0.30
b = 0.30*1.8
x = 1. + 0.011 + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 2, color = "cyan")

a = 0.25
b = 0.25*1.6
x = 1. + 0.061 + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 2, color = "magenta")

a = 0.2
b = 0.2*1.4
x = 1. + 0.11 + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 2, color = "blue")

a = 0.15
b = 0.15*1.2
x = 1. + 0.161 + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 2, color = "green")

a = 0.1
b = 0.1
x = 1. + 0.211 + a * np.cos(theta)
y = b * np.sin(theta)

plt.plot (x, y, linewidth = 2, color = "red")

plt.axhline (0., color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (1., color = 'black', linewidth = 1., linestyle = 'dotted')

plt.axis ('equal')

plt.xlabel (r'$R/R_0$', fontsize = font)
plt.ylabel (r'$Z/R_0$', fontsize = font)

#plt.show ()
plt.savefig ("Figure1.pdf")
