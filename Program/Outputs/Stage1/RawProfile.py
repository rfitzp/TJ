# -*-Python-*-
# Created by fitzpatrickr on 15 Aug 2021

# Script plots Stage1 raw plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn  = 'Stage1.nc'
ds  = nc.Dataset(fn)

psin   = ds['PSI_N']
p      = ds['P']
pp     = ds['Pp']
t      = ds['T']
tp     = ds['TTp']
q      = ds['Q']

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: Stage1 Raw Profiles")

plt.subplot(2, 3, 1)
plt.xlim((0.0, 1.0))
plt.plot(psin, p)
plt.axhline(0., color='red', linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$P$')

plt.subplot(2, 3, 2)
plt.xlim((0.0, 1.0))
plt.plot(psin, t)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$T$")

plt.subplot(2, 3, 3)
plt.xlim((0.0, 1.0))
plt.plot(psin, q)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$q$")

plt.subplot(2, 3, 4)
plt.xlim((0.0, 1.0))
plt.plot(psin, pp)
plt.axhline(0., color='red', linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$P'$")

plt.subplot(2, 3, 5)
plt.xlim((0.0, 1.0))
plt.plot(psin, tp)
plt.axhline(0., color='red', linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$TT'$")

plt.tight_layout(pad=0.5)

plt.show()
