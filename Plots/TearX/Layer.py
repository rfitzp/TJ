# Layer.py

# Plots layer quantities versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/TearX/TearX.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
tauR = ds['tauR']
tauA = ds['tauA']
tauE = ds['tauE']
taup = ds['taup']
S    = ds['S']
cb   = ds['c_beta']
db   = ds['d_beta']

rres = ds['rres']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: Layer Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (0., 1.)
 
plt.plot    (r, tauR, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',      fontsize = "15")
plt.ylabel (r'$\tau_R$', fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (0., 1.)
 
plt.plot    (r, tauA, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',      fontsize = "15")
plt.ylabel (r"$\tau_A$", fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 1.)
 
plt.plot (r, tauE, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$\tau_E$")
plt.plot (r, taup, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$\tau_\phi$")
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 1.)
 
plt.plot (r, S,  color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$S$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (0., 1.)
 
plt.plot (r, cb,  color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',       fontsize = "15")
plt.ylabel (r'$c_\beta$', fontsize = "15")

plt.subplot (3, 2, 6)

plt.xlim (0., 1.)
 
plt.plot (r, db,  color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',       fontsize = "15")
plt.ylabel (r'$d_\beta$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
