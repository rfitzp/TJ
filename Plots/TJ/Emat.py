# Emat.py

# Visualizes tearing stability matrix E_k^k'
# Also shows anti-Hermitian component of E_k^k'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Emat_r']
avaci = ds['Emat_i']
bvacr = ds['Eant_r']
bvaci = ds['Eant_i']

fig = plt.figure (figsize = (10.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Tearing Stability Matrix')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (2, 2, 1)

arrax = plt.matshow (avacr, fignum = 0, cmap= ReBu)
plt.colorbar (arrax)
plt.title (r"$Re(E_{kk'})$")

plt.subplot (2, 2, 2)

ariax = plt.matshow (avaci, fignum = 0, cmap= ReBu)
plt.colorbar (ariax)
plt.title (r"$Im(E_{kk'})$")

plt.subplot (2, 2, 3)

arrax = plt.matshow (bvacr, fignum = 0, cmap= ReBu)
plt.colorbar (arrax)
plt.title (r"$Re(E_{anti\, kk'})$")

plt.subplot (2, 2, 4)

ariax = plt.matshow (bvaci, fignum = 0, cmap= ReBu)
plt.colorbar (ariax)
plt.title (r"$Im(E_{anti\,kk'})$")

plt.tight_layout()

plt.show()
