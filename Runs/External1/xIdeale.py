# xIdeale.py

# Plots poloidal harmonics of psi, Z, and Xi components of ideal eigenfunction versus r.
# User prompted for eigenfunction number.
# Locations of rational surfaces are shown.
# Plots xi instead of Xi

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Psi_e_r']
Psi_i = ds['Psi_e_i']
Z_r   = ds['Z_e_r']
Z_i   = ds['Z_e_i']
Xi_r  = ds['xi_e_r']
Xi_i  = ds['xi_e_i']
rres  = ds['r_res']
mres  = ds['m_res']

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Ideal Eigenfunction')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = Psi_r.shape[0]

j = 1

plt.subplot(2, 1, 1)

plt.xlim(0., 1.)

plt.plot (r, Psi_r[11,j,:], color = 'black',   linewidth = 2, linestyle = 'solid', label = '$m=1$')
plt.plot (r, Psi_r[12,j,:], color = 'red',     linewidth = 2, linestyle = 'solid', label = '$m=2$')
plt.plot (r, Psi_r[13,j,:], color = 'green',   linewidth = 1, linestyle = 'solid', label = '$m=3$')
plt.plot (r, Psi_r[14,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid', label = '$m=4$')
plt.plot (r, Psi_r[15,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid', label = '$m=5$')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\psi_m$)", fontsize = "15", labelpad = -7)
plt.legend (fontsize = '13')

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)

plt.plot (r[2:], Xi_r[11,j,2:], color = 'black',   linewidth = 2, linestyle = 'solid', label = '$m=1$')
plt.plot (r, Xi_r[12,j,:], color = 'red',     linewidth = 2, linestyle = 'solid', label = '$m=2$')
plt.plot (r, Xi_r[13,j,:], color = 'green',   linewidth = 1, linestyle = 'solid', label = '$m=3$')
plt.plot (r, Xi_r[14,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid', label = '$m=4$')
plt.plot (r, Xi_r[15,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid', label = '$m=5$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($\xi_m$)", fontsize = "15")
plt.legend (fontsize = '13')

plt.tight_layout ()

plt.show ()    
