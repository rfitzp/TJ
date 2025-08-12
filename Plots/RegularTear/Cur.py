# Cur.py

# Plots tearing current eigenfunction versus radius in vicinity of rational surface
# User promted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = '../../Outputs/RegularTear/RegularTear.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
psi = ds['j']
psr = np.asarray(ds['j_r'])
rs  = ds['rres']
sg  = ds['sigma']

sig = sg[0]

font = 20

nres = len(rs)
print ("rational surface = (%d .. %d)" % (1, nres))
m   = input ("rational surface number ? ")
j   = int(m) - 1

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'RegularTear Code: Tearing Current Eigenfunction')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (1, 1, 1)

plt.xlim (rs[j]-5.*sig, rs[j]+5.*sig)
plt.ylim (1.2*min(psr[:,j]), 1.2*max(psr[:,j]))
 
plt.plot    (r, psi[:,j], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r, psr[:,j], color = 'red',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,          color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs[j],       color = 'black', linewidth = 2.0, linestyle = 'dashed')

plt.axvline (rs[j]+3.*sig, color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axvline (rs[j]-3.*sig, color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = font)
plt.ylabel (r'$j$', fontsize = font)

plt.tight_layout ()

plt.show ()    
