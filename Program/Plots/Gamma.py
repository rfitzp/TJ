# gamma.py

# Plots coefficients in expansions of Psi_x and Psi_rmp at boundary in terms of ideal eigenfunctions

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator

fn  = 'TJ.nc'
ds  = nc.Dataset(fn)
gxr = ds['gammax_r']
gxi = ds['gammax_i']
gr  = ds['gamma_r']
gi  = ds['gamma_i']

J = len(gxr)

jj = np.linspace (0, J, J, endpoint = False)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Expansion of Psi_x and Psi_rmp')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.xlim (-1, J)

plt.plot    (jj, gxr, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Re(\psi_x)$')
plt.plot    (jj, gr,  color = 'green', linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Re(\psi^{rmp})$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

ax = fig.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel (r'$j$',        fontsize = "15")
plt.ylabel (r'$\gamma_r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (-1, J)

plt.plot    (jj, gxi, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Im(\psi_x)$')
plt.plot    (jj, gi,  color = 'green', linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Im(\psi^{rmp})$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

ax = fig.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel (r'$j$',        fontsize = "15")
plt.ylabel (r'$\gamma$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
