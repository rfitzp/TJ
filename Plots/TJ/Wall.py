# Wall.py

# Visualizes wall matrices, I_m^m, J_m^m', and K_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
pvacr = ds['Iw']
pvaci = ds['Jw']
rvacr = ds['Kw']

#ar = np.power(np.abs(np.asarray(pvacr)), 0.25)
#ai = np.power(np.abs(np.asarray(pvaci)), 0.25)
#br = np.power(np.abs(np.asarray(rvacr)), 0.25)

ar = np.asarray(pvacr)
ai = np.asarray(pvaci)
br = np.asarray(rvacr)

"""
J = ar.shape[0]
for j in range(J):
    for jp in range(J):
        if ar[j,jp] > 0.:
            ar[j,jp] =   pow( ar[j,jp], 0.25)
        else:
            ar[j,jp] = - pow(-ar[j,jp], 0.25)
        if ai[j,jp] > 0.:
            ai[j,jp] =   pow( ai[j,jp], 0.25)
        else:
            ai[j,jp] = - pow(-ai[j,jp], 0.25)
        if br[j,jp] > 0.:
            br[j,jp] =   pow( br[j,jp], 0.25)
        else:
            br[j,jp] = - pow(-br[j,jp], 0.25)
"""

arp = np.amax(ar)
arm = np.amin(ar)
aip = np.amax(ai)
aim = np.amin(ai)
brp = np.amax(br)
brm = np.amin(br)

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

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Wall Matrices I, J, and K')
plt.rc ('xtick', labelsize=10) 
plt.rc ('ytick', labelsize=10)

plt.subplot (2, 2, 1)

prrax = plt.matshow (ar, fignum = 0, cmap = ReBu, vmin = -armax, vmax = armax)
plt.colorbar (prrax)
plt.title (r"$I_m^{m'}$")

plt.subplot (2, 2, 2)

priax = plt.matshow (ai, fignum = 0, cmap = ReBu, vmin = -aimax, vmax = aimax)
plt.colorbar (priax)
plt.title (r"$J_m^{m'}$")

plt.subplot (2, 2, 3)

rrrax = plt.matshow (br, fignum = 0, cmap = ReBu, vmin = -brmax, vmax = brmax)
plt.colorbar (rrrax)
plt.title (r"$K_m^{m'}$")

plt.tight_layout ()

plt.show ()
