# F72.py

# Plots function F_72 versus z

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/ECE/ECE.nc'
ds   = nc.Dataset(fn)
z    = ds['z']
F72r = ds['F72_r']
F72i = ds['F72_i']

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: ECE Absorption Function')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.xlim (z[0], z[-1])
 
plt.plot    (z, F72r, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = r'$Re(F_{7/2})$')
plt.plot    (z, F72i, color = 'red',   linewidth = 2,  linestyle = 'solid', label = r'$Im(F_{7/2})$')
plt.axhline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$z$',fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()

plt.show () 
