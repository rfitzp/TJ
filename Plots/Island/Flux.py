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
f3 = np.asarray(ds['Flux_3'])
f4 = np.asarray(ds['Flux_4'])
f5 = np.asarray(ds['Flux_5'])
f6 = np.asarray(ds['Flux_6'])
f7 = np.asarray(ds['Flux_7'])

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'Island Code: Flux Surface Averages')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 2, 1)

plt.xlim (k[0], k[-1])

plt.plot (k, f0, color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$<1>$")
plt.plot (k, f1, color = 'blue',    linewidth = 2, linestyle = 'solid', label = r"$<\cos\zeta>$")
plt.plot (k, f2, color = 'red',     linewidth = 2, linestyle = 'solid', label = r"$k <\cos\zeta>/<1>$")
plt.plot (k, f3, color = 'yellow',  linewidth = 2, linestyle = 'solid', label = r"$<cos\xi>$")
plt.plot (k, f4, color = 'cyan',    linewidth = 2, linestyle = 'solid', label = r"$<sin\xi\,\sin\zeta>$")
plt.plot (k, f5, color = 'magenta', linewidth = 2, linestyle = 'solid', label = r"$<Y^2>$")
plt.plot (k, f6, color = 'brown',   linewidth = 2, linestyle = 'solid', label = r"$G_1$")

plt.axhline ( 0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (1, 2, 2)

plt.xlim (k[0], k[-1])

plt.plot (k, f7, color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$G_2$")

plt.axhline ( 0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
