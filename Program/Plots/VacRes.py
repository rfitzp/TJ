# VacRes.py

# Visualizes vacuum solution residual matrix A_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Avac_r']
avaci = ds['Avac_i']

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (1, 2, 1)

arrax = plt.matshow (avacr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(A^{mm})$")

plt.subplot (1, 2, 2)

ariax = plt.matshow (avaci, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(A^{mm})$")

plt.tight_layout()

plt.show()
