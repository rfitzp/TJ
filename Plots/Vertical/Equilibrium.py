# Equilibrium.py

# Plots safety-factor, q, magnetic shear, s, higher order shear, s_2, pressure gradient p',
# and derivate of pressure gradient, p'', versus radius, r.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Vertical/Vertical.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
q    = ds['q']
s    = ds['s']
s2   = ds['s2']
pp   = ds['pp']
ppp  = ds['ppp']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Equilibrium Profiles')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$q$',       fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, s, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$s$',       fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, s2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$s_2$",     fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r, pp, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$p_2'$",    fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (0., 1.)

plt.plot    (r, ppp, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$p_2''$",   fontsize = "15")

plt.tight_layout ()

plt.show ()    
