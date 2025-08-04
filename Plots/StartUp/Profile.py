# Profile.py

# Plots temperature and safety-factor profiles

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/StartUp/StartUp.nc'
ds = nc.Dataset (fn)
x  = ds['x']
T  = np.asarray(ds['T_e'])
q  = ds['q']
j  = (T+1.e-8)**1.5

font = 20
fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'StartUp Code: Equilibrium Profiles')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (2, 1, 1)

plt.xlim (0., 1.)
 
plt.plot    (x, T, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = "$\overline{T}$")
plt.plot    (x, j, color = 'red',  linewidth = 2,   linestyle = 'solid', label = "$\overline{j}$")
plt.axhline (0.,   color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$x$', fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)
 
plt.plot    (x, q, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (1.,   color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$x$', fontsize = font)
plt.ylabel (r'$q$', fontsize = font)

plt.tight_layout ()

plt.show () 
#plt.savefig ("Figure3.pdf")
