# Metric.py

# Plots metric functions on plasma boundary

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'TJ.nc'
ds = nc.Dataset(fn)
t  = ds['theta']
e  = ds['eta']
cm = ds['cosmu']
ce = ds['coseta']
se = ds['sineta']
rz = ds['R2grgz']
re = ds['R2grge']

tt = np.asarray(t)
tt = tt/math.pi

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 2.)

plt.plot (tt, cm, color = 'blue', linewidth = 2, linestyle = 'solid')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r'$z$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 2.)

plt.plot    (tt, ce, color = 'blue',   linewidth = 2,   linestyle = 'solid', label = r"$\cos\,\eta$")
plt.plot    (tt, se, color = 'green', linewidth = 2,   linestyle = 'solid', label = r"$\sin\,\eta$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 2.)

plt.plot (tt, rz, color = 'blue',  linewidth = 2,   linestyle = 'solid')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r'$R^2\,\hat{\nabla} \hat{r}\cdot\hat{\nabla} z$', fontsize="15")

plt.subplot (2, 2, 4)

plt.xlim (0., 2.)

plt.plot    (tt, re, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r'$R^2\,\hat{\nabla} \hat{r}\cdot\hat{\nabla} \eta$', fontsize = "15")

plt.tight_layout()

plt.show()    
