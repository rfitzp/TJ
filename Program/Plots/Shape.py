# Shape.py

# Plots values of shaping functions at plasma boundary 

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
n  = ds['n']
Hn = ds['Hna']
Vn = ds['Vna']

nn   = np.asarray(n)
nmax = nn.shape[0]
                      
fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Edge Shaping Values')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17) 

plt.subplot (1, 2, 1)

plt.xlim (0., nmax-1)

plt.plot    (n, Hn, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$n$',        fontsize = "20")
plt.ylabel (r'$H_{n\,a}$', fontsize = "20")

plt.subplot (1, 2, 2)

plt.xlim (0., nmax-1)

plt.plot    (n, Vn, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$n$',        fontsize = "20")
plt.ylabel (r'$V_{n\,a}$', fontsize = "20")
 
plt.tight_layout ()

plt.show ()    
