# Psirwm.py

# Plots resistive wall mode eigenfunction versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = '../../Outputs/Pinch/Pinch.nc'
ds     = nc.Dataset(fn)
r      = ds['r']
rv     = ds['r_v']
psip   = ds['psi_p']
psirwm = ds['psi_rwm']
ip     = ds['InputParameters']

rs = ip[14]
bw = ip[7]
dw = ip[8]

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: Resistive Wall Mode Eigenfunction')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.2*(bw+dw))

plt.plot (r,  psip,   color = 'black', linewidth = 2, linestyle = 'solid')
plt.plot (rv, psirwm, color = 'black', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (1.,    color = 'black', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (bw,    color = 'black', linewidth = 2.0, linestyle = 'solid')
plt.axvline (bw+dw, color = 'black', linewidth = 2.0, linestyle = 'solid')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$',          fontsize = fontsize)
plt.ylabel (r'$\bar{\psi}_{rwm}$', fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
