# yZSurface.py

# Plots y and Z on plasma boundary associated with no-wall ideal eigenfunction.
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'Vertical.nc'
ds    = nc.Dataset(fn)
theta = ds['theta']
Psi_r = ds['y_surface_r']
Psi_i = -np.asarray(ds['y_surface_i'])
Xi_r  = ds['Z_surface_r']
Xi_i  = -np.asarray(ds['Z_surface_i'])
w     = np.asarray(ds['delta_W'])
p     = np.asarray(ds['delta_W_p'])
v     = np.asarray(ds['delta_W_v'])

t  = np.asarray(theta);
tt = t/math.pi

j = 0

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc('xtick', labelsize = 20) 
plt.rc('ytick', labelsize = 20)

plt.subplot (1, 2, 1)

plt.xlim (0., 2.)

plt.plot (tt, Psi_i[j,:], color = 'black', linewidth = 2, linestyle = 'solid')
#plt.plot (tt, Psi_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "20")
plt.ylabel (r"$y (a.u.)$",   fontsize = "20")

plt.subplot (1, 2, 2)

plt.xlim (0., 2.)

plt.plot (tt, Xi_i[j,:], color = 'black', linewidth = 2, linestyle = 'solid')
#plt.plot (tt, Xi_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "20")
plt.ylabel (r"$Z (a.u.)$",   fontsize = "20")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure12_4.pdf")
