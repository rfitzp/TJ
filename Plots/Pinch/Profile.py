# Profile.py

# Plots equilibrium profiles versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Pinch/Pinch.nc'
ds = nc.Dataset(fn)
r  = ds['r']
s  = np.asarray(ds['sigma'])
p  = np.asarray(ds['P'])
bt = ds['B_phi']
bp = ds['B_theta']
q  = ds['q']

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: Equilibrium Profiles')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)
plt.ylim (0., 1.05*s[0])

plt.plot (r, s, color = 'black', linewidth = 2, linestyle = 'solid')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\sigma$',  fontsize = fontsize)

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)
plt.ylim (0., 1.05*p[0])

plt.plot (r, p, color = 'black', linewidth = 2, linestyle = 'solid')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$\bar{P}$', fontsize = fontsize)

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, bt, color = 'black', linewidth = 2,   linestyle = 'solid',  label = r'$\bar{B}_\varphi$')
plt.plot    (r, bp, color = 'black', linewidth = 2,   linestyle = 'dashed', label = r'$\bar{B}_\theta$')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
#plt.ylim (0., 1.05*q[0])

plt.plot (r, q, color = 'black', linewidth = 2, linestyle = 'solid')

if (q[-1] < 0.):
    plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.ylabel (r'$q$',       fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
#plt.savefig("Figure9_1.pdf")
