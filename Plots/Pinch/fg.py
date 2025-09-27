# fg.py

# Plots f and g versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Pinch/Pinch.nc'
ds = nc.Dataset(fn)
r  = ds['r']
f  = np.asarray(ds['f'])
g  = np.asarray(ds['g'])
ip = ds['InputParameters']

rs = ip[14]
qs = ip[15]

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: f and g Profiles')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 1, 1)

plt.xlim (0., 1.)
plt.ylim (0., 1.05*f[0])

plt.plot (r, f, color = 'black', linewidth = 2, linestyle = 'solid')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$f$',       fontsize = fontsize)

plt.subplot (2, 1, 2)

plt.xlim (0.,  1.)

plt.plot    (r, g, color = 'black', linewidth = 2, linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$',    fontsize = fontsize)
plt.ylabel (r'$tanh(g/10)$', fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
