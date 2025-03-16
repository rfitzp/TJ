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

fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize  =15)

nres = len(rres)

J = Psi_r.shape[0]    

plt.subplot(2, 1, 1)

j   = 0

plt.xlim(0., 1.)

plt.plot (r, Psi_r[11,j,:], color = 'black',   linewidth = 1.5, linestyle = 'solid', label = '$m=1$')
plt.plot (r, Psi_r[12,j,:], color = 'red',     linewidth = 1.5, linestyle = 'solid', label = '$m=2$')
plt.plot (r, Psi_r[13,j,:], color = 'green',   linewidth = 1.5, linestyle = 'solid', label = '$m=3$')
plt.plot (r, Psi_r[14,j,:], color = 'blue',    linewidth = 1.5, linestyle = 'solid', label = '$m=4$')
plt.plot (r, Psi_r[15,j,:], color = 'cyan',    linewidth = 1.5, linestyle = 'solid', label = '$m=5$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$r/a$',    fontsize = "15") 
plt.ylabel (r"$\psi_m$", fontsize = "15", labelpad=-0.9)
plt.legend (fontsize = '13')

plt.subplot(2, 1, 2)

j   = 1

plt.xlim(0., 1.)

plt.plot (r, Psi_r[11,j,:], color = 'black',   linewidth = 1.5, linestyle = 'solid', label = '$m=1$')
plt.plot (r, Psi_r[12,j,:], color = 'red',     linewidth = 1.5, linestyle = 'solid', label = '$m=2$')
plt.plot (r, Psi_r[13,j,:], color = 'green',   linewidth = 1.5, linestyle = 'solid', label = '$m=3$')
plt.plot (r, Psi_r[14,j,:], color = 'blue',    linewidth = 1.5, linestyle = 'solid', label = '$m=4$')
plt.plot (r, Psi_r[15,j,:], color = 'cyan',    linewidth = 1.5, linestyle = 'solid', label = '$m=5$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$r/a$',    fontsize = "15") 
plt.ylabel (r"$\psi_m$", fontsize = "15")
plt.legend (fontsize = '13')

plt.tight_layout ()

plt.show ()    
