# Shear.py

# Plots magnetic shear, amd related functions, versus r.
# Compares edge magnetic shear with required value.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'Equilibrium.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
q    = ds['q_2']
s    = ds['s']
s2   = ds['s2']
qq   = ds['qq']
qqq  = ds['qqq']
pp   = ds['pp']
ppp  = ds['ppp']
para = ds['para']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Magnetic Shear')
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

plt.plot    (r, qq, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$r\,q'$",   fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, s,    color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (para[1], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$s$',       fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r, s2, color = 'blue',   linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$s_2$",     fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (0., 1.)

plt.plot    (r, pp, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$p'$",      fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (3, 2, 6)

plt.xlim (0., 1.)

plt.plot    (r, ppp, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$p''$",     fontsize = "15")

plt.tight_layout ()

plt.show ()    
