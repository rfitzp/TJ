# H.py

# Plots spectral convolution functions of 1st harmonic O-mode and 2nd harmonic X-mode

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/ECE/ECE.nc'
ds   = nc.Dataset (fn)
In   = ds['InputParameters']
wcO  = ds['wc']
alpO = ds['H_1^O']
alpX = ds['H_2^X']

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: ECE Spectral Convolution Functions')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (2, 1, 1)

plt.xlim (wcO[0], wcO[-1])
 
plt.plot    (wcO, alpO, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,        color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (In[5],     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\omega_c/\omega_p$',       fontsize = font)
plt.ylabel (r'$\omega_p\,H_{1}^{\,(O)}$', fontsize = font)

plt.subplot (2, 1, 2)

plt.xlim (wcO[0], wcO[-1])
 
plt.plot    (wcO, alpX, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,        color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (In[5],     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\omega_c/\omega_p$',       fontsize = font)
plt.ylabel (r'$\omega_p\,H_{2}^{\,(X)}$', fontsize = font)

plt.tight_layout ()

plt.show () 
