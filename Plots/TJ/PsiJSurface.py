# PsiJSurface.py

# Plots psi and J on plasma boundary associated with ideal eigenfunction.
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = '../../Outputs/TJ/TJ.nc'
ds     = nc.Dataset(fn)
theta  = ds['theta']
Psi_r  = ds['Psi_surface_r']
Psi_i  = ds['Psi_surface_i']
J_r    = ds['J_surface_r']
J_i    = ds['J_surface_i']
PsiJ_r = ds['PsiJ_surface_r']
PsiJ_i = ds['PsiJ_surface_i']

t  = np.asarray(theta);
tt = t/math.pi

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Ideal Eigenfunction Boundary Data')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = Psi_r.shape[0]
print ("solution number = (%d .. %d)" % (0, J-1))
j   = input ("solution number ? ")

plt.subplot (3, 1, 1)

plt.xlim (0., 2.)

plt.plot (tt, Psi_r[j,:], color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, Psi_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$\psi$",       fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim (0., 2.)

plt.plot (tt, J_r[j,:], color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, J_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.ylabel (r"$J$",          fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (0., 2.)

plt.plot (tt, PsiJ_r[j,:], color = 'blue', linewidth = 2, linestyle = 'solid')
plt.plot (tt, PsiJ_i[j,:], color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$',   fontsize = "15")
plt.ylabel (r"$\psi^\ast\,J$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
