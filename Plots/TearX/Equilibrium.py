# Equilibrium.py

# Plots equilibrium quantities versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/TearX/TearX.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
q    = ds['q']
s    = ds['s']
j    = ds['J']
jp   = ds['Jp']
l    = ds['lambda']
rres = ds['rres']

rr = np.asarray(r)
r1 = 1. + rr
lr = np.log(r1)/math.log(2.)
ss = np.asarray(s)
s1 = 1. + ss
ls = np.log(s1)

mjp = - np.asarray(jp)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: Equilibrium Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)
 
plt.plot (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$q$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)
 
plt.plot    (lr, ls, color = 'blue',  linewidth = 2,   linestyle = 'solid')
for rs in rres:
    plt.axvline (math.log(1. + rs)/math.log(2.), color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$\ln(1+r)/\ln 2$', fontsize = "15")
plt.ylabel (r'$\ln(1+s)$', fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)
 
plt.plot    (r, j,   color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r'$J$')
plt.plot    (r, mjp, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$-J'$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r'$J$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
 
plt.plot    (r, l, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',       fontsize = "15")
plt.ylabel (r'$\lambda$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
