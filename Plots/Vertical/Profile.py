# Profile.py

# Plots P1, P1, Sigma, and S5

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Vertical/Vertical.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
p1   = ds['P1']
p2   = ds['P2']
s    = ds['Sig']
s5   = ds['S5']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Profile Functions')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, p1, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$P_1$',     fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, p2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$P_2$',     fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, s, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$\Sigma$",  fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r, s5, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylabel (r"$S_5$",     fontsize = "15")
plt.xlabel (r'$\hat{r}$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
