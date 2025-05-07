# dTh.py

# Plots harmonics of temperature perturbation versus x

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
x   = np.asarray(ds['X'])
dT  = np.asarray(ds['delta_T_h'])
p   = np.asarray(ds['InputParameters'])
g    = np.asarray(ds['para'])

delta = p[5]
Finf  = g[1]

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: Harmonics of Temperature Perturbation')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

nharm = dT.shape[0]

plt.subplot (1, 1, 1)

dT0 = (dT[0,:] - x)/3.
dT1 = dT[1,:] + delta/8.**0.5

ymax = 1.1 * max (dT1)
ymin = 1.1 * min (dT1)

plt.xlim (-1.5, 1.5)
plt.ylim (ymin, ymax)

plt.plot (x, dT0,     color = 'red',     linewidth = 2, linestyle = 'solid', label = r"$\nu=0$")
plt.plot (x, dT1,     color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$\nu=1$")
plt.plot (x, dT[2,:], color = 'blue',    linewidth = 2, linestyle = 'solid', label = r"$\nu=2$")
plt.plot (x, dT[3,:], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = r"$\nu=3$")
plt.plot (x, dT[4,:], color = 'cyan',    linewidth = 2, linestyle = 'solid', label = r"$\nu=4$")
plt.plot (x, dT[5,:], color = 'magenta', linewidth = 2, linestyle = 'solid', label = r"$\nu=5$")

plt.axhline ( 0.,                  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.0 + delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/W$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

plt.show () 
#plt.savefig ("dTh.pdf")
