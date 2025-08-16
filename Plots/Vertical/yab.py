# deltaW.py

# Plots no-wall delta_W versus solution number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator

fn = '../../Outputs/Vertical/Vertical.nc'
ds = nc.Dataset(fn)
yar  = ds['ya_r']
yai  = ds['ya_i']
ybr  = ds['yb_r']
ybi  = ds['yb_i']
mpol = ds['mpol']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: y_a and y_b values')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.xlim (mpol[0]-0.5, mpol[-1]+0.5)

plt.plot    (mpol, yar, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Re(y_a)$')
plt.plot    (mpol, ybr, color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Re(y_b)$')
plt.axhline (0.,        color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (mpol[0]-0.5, mpol[-1]+0.5)

plt.plot    (mpol, yai, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Im(y_a)$')
plt.plot    (mpol, ybi, color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$Im(y_b)$')
plt.axhline (0.,        color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
