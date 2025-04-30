# Delta.py

# Plots DeltaO and DeltaX versus D

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
w  = np.asarray(ds['D'])
do = np.asarray(ds['DeltaO'])
dx = np.asarray(ds['DeltaX'])
p  = ds['InputParameters']

D     = p[12]
delta = p[5]
W     = p[13]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: DeltaO and DeltaX')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (w[0], w[-1])

plt.plot (w, do, color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$\Delta_O\,W$")
plt.plot (w, dx, color = 'blue',    linewidth = 2, linestyle = 'solid', label = r"$\Delta_X\,W$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
#plt.axvline (-delta*W/8.**0.5, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$D/\sigma$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
