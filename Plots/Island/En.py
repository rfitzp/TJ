# En.py

# Plots En functions versus k

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Island/Island.nc'
ds  = nc.Dataset(fn)
k   = np.asarray(ds['k'])
En  = np.asarray(ds['E_n'])

Nb = En.shape[0]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: En(k)')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

E0 = []
for x in En[0,:]:
    E0.append (x - 1.)

plt.subplot (1, 1, 1)

plt.xlim (k[0],   k[-1])
#plt.ylim (0,  1.1*F[-1])

plt.plot (k, E0,      color = 'red',     linewidth = 2, linestyle = 'solid', label = r'$n=0$')
plt.plot (k, En[1,:], color = 'green',   linewidth = 2, linestyle = 'solid', label = r'$n=1$')
plt.plot (k, En[2,:], color = 'blue',    linewidth = 2, linestyle = 'solid', label = r'$n=2$')
plt.plot (k, En[3,:], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = r'$n=3$')
plt.plot (k, En[4,:], color = 'cyan',    linewidth = 2, linestyle = 'solid', label = r'$n=4$')
plt.plot (k, En[5,:], color = 'magenta', linewidth = 2, linestyle = 'solid', label = r'$n=5$')
plt.axhline (0.,             color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (math.pi/2.-1.,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
