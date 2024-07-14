# Psix.py

# Plots poloidal harmonics of RMP at plasma boundary

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'TJ.nc'
ds = nc.Dataset(fn)
m  = ds['mpol']
Pr = ds['Psi^x_r']
Pi = ds['Psi^x_i']
Xr = ds['Xi_r']
Xi = ds['Xi_i']
Ur = ds['Upsilon_r']
Ui = ds['Upsilon_i']

mm = np.asarray(m).astype(float)
                      
fig = plt.figure (figsize = (9.5, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: RMP Harmonics at Plasma Boundary')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.plot    (mm, Pr, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',              fontsize = "15")
plt.ylabel (r'$real(\Psi_m^x)$', fontsize = "15")

plt.subplot (3, 2, 2)

plt.plot    (mm, Pi, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',              fontsize = "15")
plt.ylabel (r'$imag(\Psi_m^x)$', fontsize = "15")

plt.subplot (3, 2, 3)

plt.plot    (mm, Xr, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',           fontsize = "15")
plt.ylabel (r'$real(\Xi_m)$', fontsize = "15")

plt.subplot (3, 2, 4)

plt.plot    (mm, Xi, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',           fontsize = "15")
plt.ylabel (r'$imag(\Xi_m)$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.plot    (mm, Ur, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',                fontsize = "15")
plt.ylabel (r'$real(\Upsilon_m)$', fontsize = "15")

plt.subplot (3, 2, 6)

plt.plot    (mm, Ui, color = 'blue',  linewidth = 1,   linestyle = 'dotted',  marker = 'o', markersize = 3)
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$m$',                fontsize = "15")
plt.ylabel (r'$imag(\Upsilon_m)$', fontsize = "15")
 
plt.tight_layout ()

plt.show ()    
