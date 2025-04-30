# F.py

# Plots F versus k

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
k  = np.asarray(ds['k'])
F  = np.asarray(ds['F'])

G = k/2 - F

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: F(k)')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (k[0],   k[-1])
plt.ylim (0,  1.1*F[-1])

plt.plot (k, F, color = 'blue', linewidth = 2, linestyle = 'solid', label = r'$F$')
plt.plot (k, G, color = 'red',  linewidth = 2, linestyle = 'solid', label = r'$2/k - F$')
plt.axhline (3.4470e-01,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
