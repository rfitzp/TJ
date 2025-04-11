# dTo.py

# Plots temperature perturbation versus x across O- and X-points

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
x   = ds['X']
dTo = ds['delta_T_o']
dTx = ds['delta_T_x']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Temperature Perturbation across O- and X-Points')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (x[0], x[-1])

plt.plot    (x, dTo, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = 'O-point')
plt.plot    (x, dTx, color = 'red',   linewidth = 2,   linestyle = 'solid', label = 'X-point')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/W$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
