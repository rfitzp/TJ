# Rmp1.py

# Plots kth poloidal harmonic of psi and Z components of ideal response to RMP versus r.
# User prompted for poloidal harmomic mode number.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Psi_rmp_r']
Psi_i = ds['Psi_rmp_i']
Z_r   = ds['Z_rmp_r']
Z_i   = ds['Z_rmp_i']
rres  = ds['r_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Ideal Response to RMP')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

print ("poloidal mode number = (%d .. %d)" % (mpol[0], mpol[-1]))
mp   = input ("m  ? ")
jp   = int(mp) - mpol[0]

J = Psi_r.shape[0]    

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot (r, Psi_r[jp,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\psi_m$)", fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim(0., 1.)

plt.plot(r, Psi_i[jp,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Im($\psi_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot (r, Z_r[jp,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($Z_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.plot (r, Z_i[jp,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($Z_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    