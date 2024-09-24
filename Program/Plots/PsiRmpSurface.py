# PsiRmpSurface.py

# Plots psi_rmp and psi_x on plasma boundary 
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = 'TJ.nc'
ds     = nc.Dataset(fn)
theta  = ds['theta']
Psi_r  = ds['Psi_rmp_surface_r']
Psi_i  = ds['Psi_rmp_surface_i']
PsiJ_r = ds['Psi_x_surface_r']
PsiJ_i = ds['Psi_x_surface_i']

t  = np.asarray(theta);
tt = t/math.pi

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: RMP Boundary Data')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (0., 2.)

plt.plot (tt, Psi_r, color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, Psi_i, color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$\psi^{rmp}$", fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 2.)

plt.plot (tt, PsiJ_r, color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, PsiJ_i, color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$\psi_x$",     fontsize = "15")

plt.tight_layout ()

plt.show ()    
