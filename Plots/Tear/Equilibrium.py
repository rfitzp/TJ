# Equilibrium.py

# Plots equilibrium quantites versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Tear/Tear.nc'
ds = nc.Dataset(fn)
r  = ds['r']
q  = ds['q']
s  = ds['s']
j  = ds['J']
jp = ds['Jp']
l  = ds['lambda']
p  = ds['CalculationParameters']
rs = p[0];

mjp = - np.asarray(jp)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEAR Code: Equilibrium Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)
 
plt.plot    (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (1.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (3.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs,   color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$q$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)
 
plt.plot    (r, s, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs,   color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$s$', fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)
 
plt.plot    (r, j,   color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r'$J$')
plt.plot    (r, mjp, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$-J'$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs,     color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$J$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
 
plt.plot    (r, l, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs,   color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$',       fontsize = "15")
plt.ylabel (r'$\lambda$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
