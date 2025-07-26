# F.py

# Plots spatial convolution functions of 1st harmonic O-mode and 2nd harmonic X-mode

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/ECE/ECE.nc'
ds   = nc.Dataset (fn)
In   = ds['InputParameters']
wcO  = ds['R']
alpO = ds['F_1^O']
alpX = ds['F_2^X']
blpO = ds['F_1^O_fit']
blpX = ds['F_2^X_fit']

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: ECE Spatial Convolution Functions')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (2, 1, 1)

plt.xlim (wcO[-1], wcO[0])
 
plt.plot    (wcO, alpO, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.plot    (wcO, blpO, color = 'red',   linewidth = 2,  linestyle = 'dotted')
plt.axhline (0.,        color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (In[4],     color = 'black', linewidth = 1., linestyle = 'dashed')
plt.axvline (In[3],     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$R$',             fontsize = font)
plt.ylabel (r'$F_{1}^{\,(O)}$', fontsize = font)

plt.subplot (2, 1, 2)

plt.xlim (wcO[-1], wcO[0])
 
plt.plot    (wcO, alpX, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.plot    (wcO, blpX, color = 'red',   linewidth = 2,  linestyle = 'dotted')
plt.axhline (0.,        color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (In[4],     color = 'black', linewidth = 1., linestyle = 'dashed')
plt.axvline (In[3],     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$R$',             fontsize = font)
plt.ylabel (r'$F_{2}^{\,(X)}$', fontsize = font)

plt.tight_layout ()

plt.show () 
