# Hmat.py

# Visualizes vacuum response matrix H_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Hmat_r']
avaci = ds['Hmat_i']

hr = np.asarray(avacr)
hi = np.asarray(avaci)

hrp = np.amax(hr)
hrm = np.amin(hr)

if (hrp > -hrm):
    hmax = hrp
else:
    hmax = -hrm

fig = plt.figure (figsize = (12.0, 5.0))
fig.canvas.manager.set_window_title (r'TJ Code: Vacuum Response Matrix')
plt.rc ('xtick', labelsize=14) 
plt.rc ('ytick', labelsize=14)

plt.subplot (1, 2, 1)

arrax = plt.matshow (hr, fignum = 0, cmap = ReBu, vmin = -hmax, vmax = hmax)
plt.colorbar (arrax)
plt.title (r"$Re(H_{jj'})$", size=14)

plt.subplot (1, 2, 2)

ariax = plt.matshow (hi, fignum = 0, cmap = ReBu)
plt.colorbar (ariax)
plt.title (r"$Im(H_{jj'})$", size=14)

plt.tight_layout ()

plt.savefig ("Figure4.pdf")

