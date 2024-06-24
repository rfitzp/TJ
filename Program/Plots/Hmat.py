# Hmat.py

# Visualizes homogeneous vacuum response matrix, H_m^m'.
# Also shows hermitian test for matrix and hermitianized version of matrix.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
hmatr = ds['Hmat_r']
hmati = ds['Hmat_i']
hresr = ds['Hres_r']
hresi = ds['Hres_i']
hsymr = ds['Hsym_r']
hsymi = ds['Hsym_i']

fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (3, 2, 1)

arrax = plt.matshow (hmatr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(H_{mm})$")

plt.subplot (3, 2, 2)

ariax = plt.matshow (hmati, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(H_{mm})$")

plt.subplot (3, 2, 3)

brrax = plt.matshow (hresr, fignum=0)
plt.colorbar (brrax)
plt.title (r"$Re(H_{mm}-H_{mm}^\dag)$")

plt.subplot (3, 2, 4)

briax = plt.matshow (hresi, fignum=0)
plt.colorbar (briax)
plt.title (r"$Im(H_{mm}-H_{mm}^\dag)$")

plt.subplot (3, 2, 5)

arrax = plt.matshow (hsymr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(H_{mm})$")

plt.subplot (3, 2, 6)

ariax = plt.matshow (hsymi, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(H_{mm})$")

plt.tight_layout()

plt.show()
