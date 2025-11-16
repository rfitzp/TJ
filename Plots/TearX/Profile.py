# Profile.py

# Plots density and temperature profiles versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/TearX/TearX.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
ne   = ds['ne']
nep  = ds['dnedr']
Te   = ds['Te']
Tep  = ds['dTedr']
Ti   = ds['Ti']
Tip  = ds['dTidr']

rres = ds['rres']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: Density and Temperature Profiles')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)
 
plt.plot (r, ne, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',   fontsize = "15")
plt.ylabel (r'$n_e$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)
 
plt.plot (r, Te, color = 'red',  linewidth = 2,    linestyle = 'solid', label = "$T_e$")
plt.plot (r, Ti, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$T_i$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)
 
plt.plot    (r, nep, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',    fontsize = "15")
plt.ylabel (r"$n_e'$", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
 
plt.plot    (r, Tep, color = 'red',  linewidth = 2,   linestyle = 'solid', label = "$T_e'$")
plt.plot    (r, Tip, color = 'blue', linewidth = 2,   linestyle = 'solid', label = "$T_i'$")
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
