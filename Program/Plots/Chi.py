# Chi.py

# Plots RMP drives at rational surfaces

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'TJ.nc'
ds = nc.Dataset(fn)
Cr = np.asarray(ds['Chi_r'])
Ci = np.asarray(ds['Chi_i'])

nres = Cr.shape[0]
kres = np.arange (1, nres+1, 1)

A = np.zeros(nres)
T = np.zeros(nres)
for i in range (nres):
    A[i] = (Cr[i]*Cr[i] + Ci[i]*Ci[i])**0.5
    T[i] = math.atan2 (Ci[i], Cr[i]) /math.pi

fig = plt.figure (figsize = (8.5, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: RMP Drive at Rational Surfaces')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (2, 2, 1)

plt.plot    (kres, Cr, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.xlabel (r'$k$',            fontsize = "20")
plt.ylabel (r'$real(\chi_k)$', fontsize = "20")

plt.subplot (2, 2, 2)

plt.plot    (kres, Ci, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.xlabel (r'$k$',            fontsize = "20")
plt.ylabel (r'$imag(\chi_k)$', fontsize = "20")

plt.subplot (2, 2, 3)

plt.plot    (kres, A, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.xlabel (r'$k$',        fontsize = "20")
plt.ylabel (r'$|\chi_k|$', fontsize = "20")

plt.subplot (2, 2, 4)

plt.ylim (-1., 1.)

plt.plot    (kres, T, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.xlabel (r'$k$',              fontsize = "20")
plt.ylabel (r'$arg(\chi_k)/\pi$', fontsize = "20")

plt.tight_layout ()

plt.show ()    
