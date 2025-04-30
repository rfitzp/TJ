# dTh.py

# Plots particular harmonic of temperature perturbation versus x.
# User prompted for harmonic.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
x   = np.asarray(ds['X'])
dT  = np.asarray(ds['delta_T_h'])
p   = np.asarray(ds['InputParameters'])

delta = p[5];

nharm = dT.shape[0] - 1

m = input ("Harmonic number (0 .. %d) ? " % nharm)
k = int(m)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: Harmonic of Temperature Perturbation')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (-1.5, 1.5)

plt.plot (x, dT[k,:], color = 'blue',    linewidth = 2, linestyle = 'solid')

plt.axhline ( 0.,                  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.0 + delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/W$',          fontsize = "15")
plt.ylabel (r'$\delta T_\nu$', fontsize = "15")

plt.tight_layout ()

plt.show () 
