# Equilibrium.py

# Plots equilibrium quantites versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Kink/Kink.nc'
ds   = nc.Dataset(fn)
p    = ds['para']
r    = ds['r']
q    = ds['q']
s    = ds['sigma']

rs = p[0];

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'KINK Code: Equilibrium Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.xlim (0., 1.)
 
plt.plot    (r, q, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axvline (rs, color   = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$q$', fontsize = "15")

plt.tight_layout ()

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)
 
plt.plot    (r, s, color = 'red',   linewidth = 2, linestyle = 'solid')
plt.axvline (rs, color   = 'black', linewidth = 2, linestyle = 'dotted')
plt.axhline (0., color   = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$r$',      fontsize = "15")
plt.ylabel (r'$\sigma$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
