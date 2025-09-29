# FGH.py

# Plots F. G, and H versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Pinch/Pinch.nc'
ds = nc.Dataset(fn)
r  = ds['r']
f  = np.asarray(ds['F'])
g  = np.asarray(ds['G'])
h  = np.asarray(ds['H'])
ip = ds['InputParameters']

rs = ip[14]
qs = ip[15]

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: f and g Profiles')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (3, 1, 1)

plt.xlim (0., 1.)

plt.plot    (r, f, color = 'black', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$F$',       fontsize = fontsize)

plt.subplot (3, 1, 2)

plt.xlim (0.,  1.)

plt.plot    (r, g, color = 'black', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$',    fontsize = fontsize)
plt.ylabel (r'$G$', fontsize = fontsize)

plt.subplot (3, 1, 3)

plt.xlim (0.,  1.)

plt.plot    (r, h, color = 'black', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$',    fontsize = fontsize)
plt.ylabel (r'$H$', fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
