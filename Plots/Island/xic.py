# xix.py

# Plots xi_max versus X

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
x  = np.asarray(ds['X'])
c  = np.asarray(ds['xi_c'])
p  = np.asarray(ds['InputParameters'])

delta = p[5]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: \xi_c(X)')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (x[0], x[-1])
plt.ylim (0.,   1.05)

plt.plot    (x, c,                 color = 'blue', linewidth = 2,   linestyle = 'solid')
plt.axhline (1.,                   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 0.0 + delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.5 - delta/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/W$',   fontsize = "15")
plt.ylabel (r'$\xi_c$', fontsize = "15")

plt.tight_layout ()

plt.show () 
