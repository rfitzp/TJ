# Profile.py

# Plots shaping function, S_1, and profiles function, P_1, P_2, P_3, versus radius, r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn1 = '../../Outputs/Equilibrium/Equilibrium.nc'
ds1 = nc.Dataset(fn1)
r1  = ds1['r']
S1  = ds1['S1']
S2  = ds1['S2']
S3  = ds1['S3']
S4  = ds1['S4']

fn   = '../../Outputs/TJ/TJ.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
P1   = ds['P1']
P2   = ds['P2']
P3   = ds['P3']
rres = ds['r_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Shaping and Profile Functions')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r1, S1, color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = '$S_1$')
plt.plot    (r1, S2, color = 'green',  linewidth = 2,   linestyle = 'solid',  label = '$S_2$')
plt.plot    (r1, S3, color = 'yellow', linewidth = 2,   linestyle = 'solid',  label = '$S_3$')
plt.plot    (r1, S4, color = 'cyan',   linewidth = 2,   linestyle = 'solid',  label = '$S_4$')
plt.axhline (0.,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$S_n$',     fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, P1, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize = "15")
plt.ylabel(r'$P_1$',     fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, P2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$P_2$',     fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r, P3, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$P_3$',     fontsize = "15")

plt.tight_layout ()

plt.show ()    
