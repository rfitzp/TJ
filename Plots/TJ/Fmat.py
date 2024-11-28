# Fmat.py

# Visualizes resistive wall matrix F_m^m'
# Also shows anti-Hermitian component of F_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
avacr = ds['Fmtr_r']
avaci = ds['Fmtr_i']
bvacr = ds['Fmta_r']
bvaci = ds['Fmta_i']

ar = np.asarray(avacr)
ai = np.asarray(avaci)
br = np.asarray(bvacr)
bi = np.asarray(bvaci)

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
fig.canvas.manager.set_window_title (r'TJ Code: Resistive Wall Matrix F')
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12)

plt.subplot (2, 2, 1)

arrax = plt.matshow (ar, fignum = 0, cmap = ReBu, vmin = -armax, vmax = armax)
plt.colorbar (arrax)
plt.title (r"$Re(F^{mm'})$")

plt.subplot (2, 2, 2)

ariax = plt.matshow (avaci, fignum = 0, cmap = ReBu, vmin = -aimax, vmax = aimax)
plt.colorbar (ariax)
plt.title (r"$Im(F^{mm'})$")

plt.subplot (2, 2, 3)

arrax = plt.matshow (bvacr, fignum = 0, cmap = ReBu, vmin = -brmax, vmax = brmax)
plt.colorbar (arrax)
plt.title (r"$Re(F^{mm'}_{anti})$")

plt.subplot (2, 2, 4)

ariax = plt.matshow (bvaci, fignum = 0, cmap = ReBu, vmin = -bimax, vmax = bimax)
plt.colorbar (ariax)
plt.title (r"$Im(F_{anti}^{mm'})$")

plt.tight_layout ()

plt.show ()

