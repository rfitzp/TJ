# Delta.py

# Plots DeltaO and DeltaX versus D

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
w   = np.asarray(ds['D'])
do  = np.asarray(ds['DeltaO'])
dx  = np.asarray(ds['DeltaX'])
dox = (do + dx) /2.
dxx = do + dx

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: DeltaO and DeltaX')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (1, 1, 1)

plt.xlim (w[0], w[-1])

plt.plot (w, do,  color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"$G_{eccd}(0)$")
plt.plot (w, dxx, color = 'red',   linewidth = 2, linestyle = 'dotted', label = r"$2\,G_{eccd}(\pi/2)$")
plt.plot (w, dox, color = 'green', linewidth = 2, linestyle = 'solid',  label = r"$G_{eccd}(\pi/2)$")
plt.plot (w, dx,  color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"$G_{eccd}(\pi)$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$d/D$', fontsize = fontsize)
plt.legend (fontsize = fontsize);

plt.tight_layout ()

#plt.show () 
plt.savefig("W40.pdf")
