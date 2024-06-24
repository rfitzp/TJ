# Chi.py

# Plots resonant magnetic perturbation response vector versus poloidal mode number associated with given rational surface
# User prompted for rational surface number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
chi_r = ds['Chi_r']
chi_i = ds['Chi_i']
chi_m = ds['Chi_m']
chi_a = ds['Chi_a']
rres  = ds['rres']

fig = plt.figure (figsize=(12.0, 8.0))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15)

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m   = input ("rational surface number ? ")
j   = int(m) - 1

plt.subplot(2, 2, 1)

plt.plot(mpol, chi_r[j,:], color='blue', marker = 'o', mfc = 'r', linewidth = 1, linestyle = 'dashed')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel(r'$m$', fontsize="15")
plt.ylabel(r"real($\chi_m$)", fontsize="15")

plt.subplot(2, 2, 2)

plt.plot(mpol, chi_i[j,:], color='blue', marker = 'o', mfc = 'r', linewidth = 1, linestyle = 'dashed')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel(r'$m$', fontsize="15")
plt.ylabel(r"imag($\chi_m$)", fontsize="15")

plt.subplot(2, 2, 3)

plt.plot(mpol, chi_m[j,:], color='blue', marker = 'o', mfc = 'r', linewidth = 1, linestyle = 'dashed')

plt.xlabel(r'$m$', fontsize="15")
plt.ylabel(r"$|\chi_m|$", fontsize="15")

plt.subplot(2, 2, 4)

plt.ylim(-1.,1.);

plt.plot(mpol, chi_a[j,:], color='blue', marker = 'o', mfc = 'r', linewidth = 1, linestyle = 'dashed')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel(r'$m$', fontsize="15")
plt.ylabel(r"arg($\chi_m)/\pi$", fontsize="15")

plt.tight_layout()

plt.show()    
