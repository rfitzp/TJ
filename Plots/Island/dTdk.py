# dTdk.py

# Plots dTdk versus k

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/Island/Island.nc'
ds   = nc.Dataset(fn)
k    = np.asarray(ds['k'])
dTdk = np.asarray(ds['dTdk'])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: dTdk(k)')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (k[0], k[-1])

plt.plot    (k, dTdk, color = 'blue', linewidth = 2, linestyle = 'solid')
plt.axhline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$',    fontsize = "15")
plt.ylabel (r'$dTdk$', fontsize = "15")

plt.tight_layout ()

plt.show () 
