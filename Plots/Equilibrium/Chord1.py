# Chord1.py

# Plots metric quantities along tilted central chord

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn1 = '../../Outputs/Equilibrium/Equilibrium.nc'
ds1 = nc.Dataset(fn1)
dRdr = ds1['dRdr_eq']
dRdt = ds1['dRdt_eq']
dZdr = ds1['dZdr_eq']
dZdt = ds1['dZdt_eq']

Nf = len(dRdr)

x = np.linspace (0., 1., Nf);

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Metric Quantities Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot    (x, dRdr, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',                          fontsize = "15")
plt.ylabel (r'$\partial R/\partial\hat{r}$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

plt.plot    (x, dRdt, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',                                        fontsize = "15")
plt.ylabel (r'$(\partial R/\partial\hat{\theta})/\hat{r}$', fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot    (x, dZdr, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',                          fontsize = "15")
plt.ylabel (r'$\partial Z/\partial\hat{r}$', fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.plot    (x, dZdt, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',                                        fontsize = "15")
plt.ylabel (r'$(\partial Z/\partial\hat{\theta})/\hat{r}$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
