# P.py

# Plots pressure and pressure gradient profiles versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Pinch/Pinch.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
s   = np.asarray(ds['sigma'])
p   = np.asarray(ds['P'])
pc  = np.asarray(ds['P_c'])
pp  = -np.asarray(ds['dPdr'])
ppc = -np.asarray(ds['dPdr_crit'])

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: Pressure Profiles')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 1, 2)

plt.xlim (0., 0.2)
plt.ylim (0.023, 0.026)

plt.plot (r, pc, color = 'black', linewidth = 2, linestyle = 'solid')
plt.plot (r, p,  color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\bar{P}$', fontsize = fontsize)

plt.subplot (2, 1, 1)

pmax = max (pp)

plt.xlim (0., 1.)
plt.ylim (0., 1.05*pmax)

plt.plot    (r, pp,  color = 'black', linewidth = 2,   linestyle = 'solid')
plt.plot    (r, ppc, color = 'black', linewidth = 2,   linestyle = 'dashed')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$',             fontsize = fontsize)
plt.ylabel (r'$-d\bar{P}/d\bar{r}$',  fontsize = fontsize)
plt.tight_layout ()

plt.show ()    
#plt.savefig("Figure9_4.pdf")
