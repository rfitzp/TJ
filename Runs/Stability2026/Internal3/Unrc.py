# Unrc.py

# Plots poloidal harmonics of psi and Z components of unreconnected solution vector associated with given rational surface versus r.
# User prompted for rational surface number.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Psi_unrc_r']
Psi_i = ds['Psi_unrc_i']
Z_r   = ds['Z_unrc_r']
Z_i   = ds['Z_unrc_i']
T_u   = ds['Torque_unrc']
rres  = ds['r_res']
mres  = ds['m_res']

fig = plt.figure (figsize = (8.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Unreconnected Eigenfunctions')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot(2, 1, 1)

plt.xlim(0., 1.)

plt.plot (r, Psi_r[11,0,:], color = 'black',   linewidth = 2, linestyle = 'solid',  label = r'$m=1$')
plt.plot (r, Psi_r[12,0,:], color = 'black',   linewidth = 2, linestyle = 'dashed', label = r'$m=2$')
plt.plot (r, Psi_r[13,0,:], color = 'black',   linewidth = 2, linestyle = 'dotted', label = r'$m=3$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15") 
plt.ylabel (r"$\psi_m$",  fontsize = "15")
plt.legend (fontsize = 15)

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)

plt.plot (r, Psi_r[11,1,:], color = 'black',   linewidth = 2, linestyle = 'solid',  label = r'$m=1$')
plt.plot (r, Psi_r[12,1,:], color = 'black',   linewidth = 2, linestyle = 'dashed', label = r'$m=2$')
plt.plot (r, Psi_r[13,1,:], color = 'black',   linewidth = 2, linestyle = 'dotted', label = r'$m=3$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$\psi_m$",  fontsize = "15")
plt.legend (fontsize = 15, loc = 'lower center')

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure13_10.pdf")
