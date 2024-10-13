# Xi.py

# Plots Fourier harmonics of Xu on plasma boundary associated with ideal eigenfunction.
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator

fn   = '../../Outputs/TJ/TJ.nc'
ds   = nc.Dataset(fn)
wr   = ds['Xi_e_r']
wi   = ds['Xi_e_i']
mpol = ds['mpol']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Fourier Harmonics of Xi on Plasma Boundary')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = mpol.shape[0]
print ("solution number = (%d .. %d)" % (0, J))
j   = input ("solution number ? ")

xir = wr[:,j,-1]
xii = wi[:,j,-1]

plt.subplot (2, 1, 1)

plt.ylim (-1.05, 1.05)

plt.plot    (mpol, xir,  color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 7)
plt.axhline (0.,         color = 'black', linewidth = 1.5, linestyle = 'dotted')

ax = fig.gca ()
ax.xaxis.set_major_locator (MaxNLocator(integer=True))

plt.xlabel (r'$m$',         fontsize = "15")
plt.ylabel (r'$Re(\Xi_m)$', fontsize = "15")

plt.subplot (2, 1, 2)

plt.ylim (-1.05, 1.05)

plt.plot    (mpol, xii,  color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 7)
plt.axhline (0.,         color = 'black', linewidth = 1.5, linestyle = 'dotted')

ax = fig.gca ()
ax.xaxis.set_major_locator (MaxNLocator(integer=True))

plt.xlabel (r'$m$',         fontsize = "15")
plt.ylabel (r'$Im(\Xi_m)$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
