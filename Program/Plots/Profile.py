# Profile.py

# Plots shaping function, S_1, and profiles function, P_1, P_2, P_3, versus radius, r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'TJ.nc'
ds = nc.Dataset(fn)
r  = ds['r']
S1 = ds['S1']
P1 = ds['P1']
P2 = ds['P2']
P3 = ds['P3']
rres = ds['rres']

fig = plt.figure(figsize=(12.0, 8.0))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 

plt.subplot(2, 2, 1)

plt.xlim(0., 1.)

plt.plot(r, S1, color='blue', linewidth = 2, linestyle = 'solid')
plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r'$S_1$', fontsize="15")

plt.subplot(2, 2, 2)

plt.xlim(0., 1.)

plt.plot(r, P1, color='blue', linewidth = 2, linestyle = 'solid')
plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r'$P_1$', fontsize="15")

plt.subplot(2, 2, 3)

plt.xlim(0., 1.)

plt.plot(r, P2, color='blue', linewidth = 2, linestyle = 'solid')
plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r'$P_2$', fontsize="15")

plt.subplot(2, 2, 4)

plt.xlim(0., 1.)

plt.plot(r, P3, color='blue', linewidth = 2, linestyle = 'solid')
plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r'$P_3$', fontsize="15")

plt.tight_layout()

plt.show()    
