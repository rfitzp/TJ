# Cmat.py

# Visualizes wall matrix C_m^m'
# Also shows anti-Hermitian component of C_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Cmat_r']
avaci = ds['Cmat_i']
bvacr = ds['Cant_r']
bvaci = ds['Cant_i']
mpol  = np.asarray(ds['mpol'])

ar = np.asarray(avacr)
ai = np.asarray(avaci)
br = np.asarray(bvacr)
bi = np.asarray(bvaci)
J = avacr.shape[0]
for i in range(J):
    if mpol[i] == 0:
        ar[i,i] = ar[i,i]/100.

arp = np.amax(ar)
arm = np.amin(ar)
aip = np.amax(ai)
aim = np.amin(ai)
brp = np.amax(br)
brm = np.amin(br)
bip = np.amax(bi)
bim = np.amin(bi)

if (arp > -arm):
    armax = arp
else:
    armax = -arm
if (aip > -aim):
    aimax = aip
else:
    aimax = -aim
if (brp > -brm):
    brmax = brp
else:
    brmax = -brm
if (bip > -bim):
    bimax = bip
else:
    bimax = -bim

fig = plt.figure (figsize = (10.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Ideal Wall Matrix C')
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12)

plt.subplot (2, 2, 1)

arrax = plt.matshow (ar, fignum = 0, cmap = ReBu, vmin = -armax, vmax = armax)
plt.colorbar (arrax)
plt.title (r"$Re(C^{mm'})$")

plt.subplot (2, 2, 2)

ariax = plt.matshow (avaci, fignum = 0, cmap = ReBu, vmin = -aimax, vmax = aimax)
plt.colorbar (ariax)
plt.title (r"$Im(C^{mm'})$")

plt.subplot (2, 2, 3)

arrax = plt.matshow (bvacr, fignum = 0, cmap = ReBu, vmin = -brmax, vmax = brmax)
plt.colorbar (arrax)
plt.title (r"$Re(C^{mm'}_{anti})$")

plt.subplot (2, 2, 4)

ariax = plt.matshow (bvaci, fignum = 0, cmap = ReBu, vmin = -bimax, vmax = bimax)
plt.colorbar (ariax)
plt.title (r"$Im(C_{anti}^{mm'})$")

plt.tight_layout ()

plt.show ()

