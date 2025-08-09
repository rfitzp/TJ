# Gmat.py

# Visualizes perfect-wall vacuum response matrix H_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
avacr = ds['Gmat_r']
avaci = ds['Gmat_i']
mpol  = np.asarray(ds['mpol'])

hr = np.asarray(avacr)
hi = np.asarray(avaci)
J = avacr.shape[0]
for i in range(J):
    if mpol[i] == 0:
        hr[i,i] = hr[i,i]/100.

hrp = np.amax(hr)
hrm = np.amin(hr)

if (hrp > -hrm):
    hmax = hrp
else:
    hmax = -hrm

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Perfect-Wall Vacuum Response Matrix')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (1, 2, 1)

arrax = plt.matshow (hr, fignum = 0, cmap = ReBu, vmin = -hmax, vmax = hmax)
plt.colorbar (arrax)
plt.title (r"$Re(G_{mm'})$")

plt.subplot (1, 2, 2)

ariax = plt.matshow (hi, fignum = 0, cmap = ReBu)
plt.colorbar (ariax)
plt.title (r"$Im(G_{mm'})$")

plt.tight_layout ()

plt.show ()
