# Solutions.py

# Plots torque and shielding curves associated with kth rational surface in plasma
# User prompted for k.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'layer.nc'
ds    = nc.Dataset(fn)
w     = ds['omega_r']
T     = ds['T_res']
X     = ds['Xi_res']

fig = plt.figure (figsize=(12.0, 8.0))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15)

nres = T.shape[0]
print ("rational surfaces = (%d .. %d)" % (1, nres))
m   = input ("k ? ")
k   = int(m) - 1

ww = np.asarray(w)/1.e3;

plt.subplot(2, 1, 1)

plt.plot(ww, T[k,:], color='blue', linewidth = 2, linestyle = 'solid')

plt.xlim(ww[0], ww[-1])

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel(r'$\omega$ (kHz)', fontsize="15")
plt.ylabel(r"$\delta T$", fontsize="15")

plt.subplot(2, 1, 2)

plt.plot(ww, X[k,:], color='blue', linewidth = 2, linestyle = 'solid')

plt.xlim(ww[0], ww[-1])

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1., color='black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel(r'$\omega$ (kHz)', fontsize="15")
plt.ylabel(r"$\Xi$", fontsize="15")

plt.tight_layout()

plt.show() 
