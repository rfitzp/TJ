# Ideale1.py

# Plots kth poloidal harmonic of y and Z components of ideal eigenfunction versus r.
# User prompted for eigenfunction number and poloidal harmomic mode number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['y_e_r']
Psi_i = ds['y_e_i']
Z_r   = ds['Z_e_r']
Z_i   = ds['Z_e_i']


fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Ideal Eigenfunction')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = Psi_r.shape[0]
print ("solution number = (%d .. %d)" % (0, J-1))
j = input ("solution number ? ")

print ("poloidal mode number = (%d .. %d)" % (mpol[0], mpol[-1]))
mp   = input ("m  ? ")
jp   = int(mp) - mpol[0]

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot (r, Psi_r[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($y_m$)", fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim(0., 1.)

plt.plot(r, Psi_i[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($y_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot (r, Z_r[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($Z_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.plot (r, Z_i[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($Z_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
