# Amat.py

# Visualizes vacuum matrix A_m^m'
# Also shows anti-Hermitian component of A_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Amat_r']
avaci = ds['Amat_i']
bvacr = ds['Aant_r']
bvaci = ds['Aant_i']

fig = plt.figure (figsize = (10.0, 8.0))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (2, 2, 1)

arrax = plt.matshow (avacr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(A^{mm'})$")

plt.subplot (2, 2, 2)

ariax = plt.matshow (avaci, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(A^{mm'})$")

plt.subplot (2, 2, 3)

arrax = plt.matshow (bvacr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(A^{mm'}_{anti})$")

plt.subplot (2, 2, 4)

ariax = plt.matshow (bvaci, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(A_{anti}^{mm'})$")

plt.tight_layout()

plt.show()
