# VacRes.py

# Visualizes vacuum solution residual matrices A_m^m', B_m^m', and C_m^m'.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Avac_r']
avaci = ds['Avac_i']
bvacr = ds['Bvac_r']
bvaci = ds['Bvac_i']
cvacr = ds['Cvac_r']
cvaci = ds['Cvac_i']

fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (3, 2, 1)

arrax = plt.matshow (avacr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(A^{mm})$")

plt.subplot (3, 2, 2)

ariax = plt.matshow (avaci, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(A^{mm})$")

plt.subplot (3, 2, 3)

brrax = plt.matshow (bvacr, fignum=0)
plt.colorbar (brrax)
plt.title (r"$Re(B^{mm})$")

plt.subplot (3, 2, 4)

briax = plt.matshow (bvaci, fignum=0)
plt.colorbar (briax)
plt.title (r"$Im(B^{mm})$")

plt.subplot (3, 2, 5)

crrax = plt.matshow (cvacr, fignum=0)
plt.colorbar (crrax)
plt.title (r"$Re(C^{mm})$")

plt.subplot (3, 2, 6)

criax = plt.matshow (cvaci, fignum=0)
plt.colorbar (criax)
plt.title (r"$Im(C^{mm})$")

plt.tight_layout()

plt.show()
