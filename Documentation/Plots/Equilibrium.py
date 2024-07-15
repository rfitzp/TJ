# Equilibrium.py

# Plots components of aspect-ratio expanded equilibrium versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
r  = ds['r']
p2 = ds['p_2']
q  = ds['q_2']

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Equilibrium Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, p2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$p_2$',     fontsize = "15")

plt.subplot (1, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, q , color = 'blue',   linewidth = 2,   linestyle = 'solid', label = r"$q_2$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (3.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (4.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r'$q$',     fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.tight_layout ()

plt.savefig ("Figure1.pdf")    
