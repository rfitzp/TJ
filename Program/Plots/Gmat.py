# Gmat.py

# Visualizes inhomogeneous vacuum response matrix, G_m^m'.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
gmatr = ds['Gmat_r']
gmati = ds['Gmat_i']

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (1, 2, 1)

arrax = plt.matshow (gmatr, fignum=0)
plt.colorbar (arrax)
plt.title (r"$Re(G_{mm})$")

plt.subplot (1, 2, 2)

ariax = plt.matshow (gmati, fignum=0)
plt.colorbar (ariax)
plt.title (r"$Im(G_{mm})$")

plt.tight_layout()

plt.show()
