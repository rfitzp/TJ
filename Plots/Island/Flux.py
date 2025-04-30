# Flux.py

# Plots flux surface averages versus k

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
k  = np.asarray(ds['kk'])
f0 = np.asarray(ds['Flux_0'])
f1 = np.asarray(ds['Flux_1'])
f2 = np.asarray(ds['Flux_2'])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: Flux Surface Averages')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (k[0], k[-1])

plt.plot (k, f0, color = 'green', linewidth = 2, linestyle = 'solid', label = r"$<1>$")
plt.plot (k, f1, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$<\cos\zeta>$")
plt.plot (k, f2, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$k <\cos\zeta>/<1>$")

plt.axhline ( 0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
