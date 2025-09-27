# Psi.py

# Plots ideal eigenfunctions Psi_p, Psi_nw, and Psi_pw versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/Pinch/Pinch.nc'
ds    = nc.Dataset(fn)
r     = ds['r']
rv    = ds['r_v']
psip  = ds['psi_p']
psinw = ds['psi_nw']
psipw = ds['psi_pw']
ip    = ds['InputParameters']

rs = ip[14]
bw = ip[7]
dw = ip[8]

fontsize = 20

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Pinch Code: Ideal Eigenfunctions')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.2*bw)

plt.plot (r,  psip,  color = 'black', linewidth = 2, linestyle = 'solid',   label = r"$\psi_p$")
plt.plot (rv, psinw, color = 'black', linewidth = 2, linestyle = 'dashed',  label = r"$\psi_{nw}$")
plt.plot (rv, psipw, color = 'black', linewidth = 2, linestyle = 'dashdot', label = r"$\psi_{pw}$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (1.,    color = 'black', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (bw,    color = 'black', linewidth = 1.5, linestyle = 'solid')
plt.axvline (bw+dw, color = 'black', linewidth = 1.5, linestyle = 'solid')

if rs > 0.:
    plt.axvline (rs, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\bar{r}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
