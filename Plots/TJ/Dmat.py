# Wmat.py

# Visualizes plasma energy matrix D

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
bvacr = ds['Dmat_r']
bvaci = ds['Dmat_i']

br = np.asarray(bvacr)
bi = np.asarray(bvaci)

brp = np.amax(br)
brm = np.amin(br)
bip = np.amax(bi)
bim = np.amin(bi)

if (brp > -brm):
    brmax = brp
else:
    brmax = -brm
if (bip > -bim):
    bimax = bip
else:
    bimax = -bim

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Plasma Energy Matrix D')
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12)

plt.subplot (1, 2, 1)

arrax = plt.matshow (br, fignum = 0, cmap = ReBu, vmin = -brmax, vmax = brmax)
plt.colorbar (arrax)
plt.title (r"$Re(D^{mm'})$")

plt.subplot (1, 2, 2)

ariax = plt.matshow (bi, fignum = 0, cmap = ReBu, vmin = -bimax, vmax = bimax)
plt.colorbar (ariax)
plt.title (r"$Im(D^{mm'})$")

plt.tight_layout ()

plt.show ()

