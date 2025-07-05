# Tunrc.py

# Plots angular momentum flux associated with pair of unreconnected solutions versus r
# User prompted for rational surface numbers, k, kp.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
T_f   = ds['Torque_pair_unrc']
rres  = ds['rres']
mres  = ds['m_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Unreconnected Electromagnetic Torques')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

nres = len(rres)
print ("rational surfaces = (%d .. %d)" % (1, nres))
m   = input ("rational surface number 1 ? ")
j   = int(m) - 1

mp   = input ("rational surface number 2 ? ")
jp   = int(mp) - 1

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)

plt.plot (r, T_f[j,jp,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "19")
plt.ylabel (r"$T_\phi$",  fontsize = "19")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure9.pdf")
