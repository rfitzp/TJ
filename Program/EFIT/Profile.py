# Profile.py

# Plots EFIT profiles versus Psi_N

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn  = 'EFIT.nc'
ds  = nc.Dataset(fn)

psin   = ds['PSI_N']
p      = ds['P']
px     = ds['Pp']
t      = ds['T']
tp     = ds['TTp']
q      = ds['Q']

pp  = np.asarray(p)/1.e6
ppx = np.asarray(px)/1.e6

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("EFIT: Equilibrium Profiles")

plt.subplot(2, 3, 1)

plt.xlim((0.0, 1.0))

plt.plot(psin, pp, linewidth = 1, color = 'blue')

plt.axhline(0., color = 'red', linewidth = 1, linestyle = 'dotted')

plt.xlabel('$\\Psi_N$', fontsize = 15)
plt.ylabel('$P$ (MPa)',   fontsize = 15)

plt.subplot(2, 3, 2)

plt.xlim((0.0, 1.0))

plt.plot(psin, t, linewidth = 1, color = 'blue')

plt.xlabel('$\\Psi_N$', fontsize = 15)
plt.ylabel("$T$ (Tm)",   fontsize = 15)

plt.subplot(2, 3, 3)

plt.xlim((0.0, 1.0))

plt.plot(psin, q, linewidth = 1, color = 'blue')

plt.xlabel('$\\Psi_N$', fontsize = 15)
plt.ylabel("$q$",       fontsize = 15)

plt.subplot(2, 3, 4)

plt.xlim((0.0, 1.0))

plt.plot(psin, ppx, linewidth = 1, color = 'blue')

plt.axhline(0., color = 'red', linewidth = 1, linestyle = 'dotted')

plt.xlabel('$\\Psi_N$', fontsize = 15)
plt.ylabel("$P' (MA/m^3)$",      fontsize = 15)

plt.subplot(2, 3, 5)

plt.xlim((0.0, 1.0))

plt.plot(psin, tp, linewidth = 1, color = 'blue')

plt.axhline(0., color = 'red', linewidth = 1, linestyle = 'dotted')

plt.xlabel('$\\Psi_N$', fontsize = 15)
plt.ylabel("$TT'$ (T)",  fontsize = 15)

plt.tight_layout()

plt.show()
