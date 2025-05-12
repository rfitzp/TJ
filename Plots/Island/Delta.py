# Delta.py

# Plots Geccd versus W

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
w   = np.asarray(ds['W'])
do  = np.asarray(ds['DeltaO'])
dx  = np.asarray(ds['DeltaX'])
dox = (do + dx)/2.

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: G_eccd')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (0., w[-1])

plt.plot (w, do,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$G_{eccd}(0)$")
plt.plot (w, dox, color = 'green', linewidth = 2, linestyle = 'solid', label = r"$G_{eccd}(\pi/2)$")
plt.plot (w, dx,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$G_{eccd}(\pi)$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$W/D$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
