# Profile2.py

# Plots additional equilibrium profiles versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Pinch/Pinch.nc'
ds = nc.Dataset(fn)
r  = ds['r']
s  = np.asarray(ds['sigma'])
sp = np.asarray(ds['sigma_p'])
b1 = ds['beta_1']
b2 = ds['beta_2']
ip = ds['InputParameters']

rs = ip[14]
qs = ip[15]

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: Additional Equilibrium Profiles')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)
plt.ylim (0., 1.05*s[0])

plt.plot (r, s, color = 'black', linewidth = 2, linestyle = 'solid')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\sigma$',  fontsize = fontsize)

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, sp, color = 'black', linewidth = 2, linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r"$\sigma'$", fontsize = fontsize)

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, b1, color = 'black', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\beta_1$', fontsize = fontsize)

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
#plt.ylim (0., 1.05*q[0])

plt.plot    (r, b2, color = 'black', linewidth = 2, linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\beta_2$', fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
#plt.savefig("Figure9_2.pdf")
