# dT0.py

# Plots zeroth harmonic of temperature perturbation versus x

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
x  = np.asarray(ds['X'])
dT = np.asarray(ds['delta_T_h'])
p  = np.asarray(ds['InputParameters'])
g  = np.asarray(ds['para'])

delta = p[5]
Finf  = g[0]

dT0 = dT[0,:] - x + Finf 

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: Zeroth Harmonic of Temperature Perturbation')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

nharm = dT.shape[0]

plt.subplot (1, 1, 1)

plt.xlim (x[0], x[-1])
#plt.ylim (-0.15, 0.15)

plt.plot (x, dT0, color = 'red', linewidth = 2, linestyle = 'solid')

plt.axhline ( 0.,                  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.0 + delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/W$', fontsize = "15")
plt.ylabel (r'$\delta T_0 - x - \delta T_{0\,\infty}$', fontsize = "15")

plt.tight_layout ()

plt.show () 
