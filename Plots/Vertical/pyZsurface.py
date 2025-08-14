# pyZSurface.py

# Plots y and Z on plasma boundary associated with perfect-wall ideal eigenfunction.
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
theta = ds['theta']
Psi_r = ds['py_surface_r']
Psi_i = ds['py_surface_i']
Xi_r  = ds['pZ_surface_r']
Xi_i  = ds['pZ_surface_i']
w     = np.asarray(ds['pdelta_W'])
p     = np.asarray(ds['pdelta_W_p'])
v     = np.asarray(ds['pdelta_W_v'])

t  = np.asarray(theta);
tt = t/math.pi

J = Psi_r.shape[0]
print ("solution number = (%d .. %d)" % (0, J-1))
j   = int(input ("solution number ? "))

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Perfect-Wall Ideal Eigenfunction Boundary Data: deltaW = %10.3e %10.3e %10.3e' % (w[j], p[j], v[j]))
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (0., 2.)

plt.plot (tt, Psi_r[j,:], color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, Psi_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$y$",          fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 2.)

plt.plot (tt, Xi_r[j,:], color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, Xi_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$Z$",          fontsize = "15")

plt.tight_layout ()

plt.show ()    
