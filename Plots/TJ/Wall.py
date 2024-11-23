# Wall.py

# Visualizes wall matrices, R_m^m' and S_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
pvacr = ds['Rwal_r']
pvaci = ds['Rwal_i']
rvacr = ds['Swal_r']
rvaci = ds['Swal_i']

ar = np.asarray(pvacr)
ai = np.asarray(pvaci)
br = np.asarray(rvacr)
bi = np.asarray(rvaci)

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

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Wall Matrices R and S')
plt.rc ('xtick', labelsize=10) 
plt.rc ('ytick', labelsize=10)

plt.subplot (2, 2, 1)

prrax = plt.matshow (pvacr, fignum = 0, cmap = ReBu, vmin = -armax, vmax = armax)
plt.colorbar (prrax)
plt.title (r"$Re(R_m^{m'})$")

plt.subplot (2, 2, 2)

priax = plt.matshow (pvaci, fignum = 0, cmap = ReBu, vmin = -aimax, vmax = aimax)
plt.colorbar (priax)
plt.title (r"$Im(R_m^{m'})$")

plt.subplot (2, 2, 3)

rrrax = plt.matshow (rvacr, fignum = 0, cmap = ReBu, vmin = -brmax, vmax = brmax)
plt.colorbar (rrrax)
plt.title (r"$Re(S_m^{m'})$")

plt.subplot (2, 2, 4)

rriax = plt.matshow (rvaci, fignum = 0, cmap = ReBu, vmin = -bimax, vmax = bimax)
plt.colorbar (rriax)
plt.title (r"$Im(S_m^{m'})$")

plt.tight_layout ()

plt.show ()
