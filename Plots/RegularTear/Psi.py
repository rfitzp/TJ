# Psi.py

# Plots tearing eigenfunction versus radius
# User promted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/RegularTear/RegularTear.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
psi = ds['psi']
psr = ds['psi_r']
rs  = ds['rres']
sg  = ds['sigma']

sig = sg[0]

font = 20

nres = len(rs)
print ("rational surface = (%d .. %d)" % (1, nres))
m   = input ("rational surface number ? ")
j   = int(m) - 1

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'RegularTear Code: Tearing Eigenfunction')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)
 
plt.plot    (r, psi[:,j], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r, psr[:,j], color = 'red',   linewidth = 2,   linestyle = 'dashed')
plt.axhline (0.,          color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs[j],       color = 'black', linewidth = 2.0, linestyle = 'dashed')

plt.axvline (rs[j]+3.*sig, color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axvline (rs[j]-3.*sig, color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$',    fontsize = font)
plt.ylabel (r'$\psi$', fontsize = font)

plt.tight_layout ()

plt.show ()    
