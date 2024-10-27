# Psi.py

# Plots tearing eigenfunction versus radius
# User promted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/Tear/Tear.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
psi = ds['psi']
rs  = ds['rres'];

nres = len(rs)
print ("rational surface = (%d .. %d)" % (1, nres))
m   = input ("rational surface number ? ")
j   = int(m) - 1

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEAR Code: Tearing Eigenfunction')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)
 
plt.plot    (r, psi[:,j], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,          color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs[j],       color = 'black', linewidth = 2.0, linestyle = 'dashed')

plt.xlabel (r'$r$',    fontsize = "15")
plt.ylabel (r'$\psi$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
