# wPSmat.py

# Visualizes wall matrix (P^dag S - R^dag Q)_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
avacr = ds['wPSmat_r']
avaci = ds['wPSmat_i']

hr = np.asarray(avacr)
hi = np.asarray(avaci)

hrp = np.amax(hr)
hrm = np.amin(hr)

if (hrp > -hrm):
    hmax = hrp
else:
    hmax = -hrm

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Wall Matrix PS = P^dag S - R^dag Q')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12)

plt.subplot (1, 2, 1)

arrax = plt.matshow (hr, fignum = 0, cmap = ReBu, vmin = -hmax, vmax = hmax)
plt.colorbar (arrax)
plt.title (r"$Re(PS_{mm'})$")

plt.subplot (1, 2, 2)

ariax = plt.matshow (hi, fignum = 0, cmap = ReBu)
plt.colorbar (ariax)
plt.title (r"$Im(PS_{mm'})$")

plt.tight_layout ()

plt.show ()
