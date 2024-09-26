# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage1 plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn  = 'Stage1.nc'
ds  = nc.Dataset(fn)

psin   = ds['PSI_N']
p      = ds['p']
pp     = ds['pp']
t      = ds['t']
tp     = ds['ttp']
q      = ds['q']

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: Stage1 Profiles")

plt.subplot(2, 3, 1)
plt.xlim((0.0, 1.0))
plt.plot(psin, p)
plt.axhline(0.,     color='red',   linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$P$')

plt.subplot(2, 3, 2)
plt.xlim((0.0, 1.0))
plt.plot(psin, t)
plt.axhline(1., color='red', linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$g$")

plt.subplot(2, 3, 3)
plt.xlim((0.0, 1.0))
plt.plot(psin, q)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$q$")

plt.subplot(2, 3, 4)
plt.xlim((0.0, 1.0))
plt.plot(psin, pp)
plt.axhline(0.,     color='red',   linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$P'$")

plt.subplot(2, 3, 5)
plt.xlim((0.0, 1.0))
plt.plot(psin, tp)
plt.axhline(0.,     color='red',   linewidth=1, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$gg'$")

plt.tight_layout(pad=0.5)

plt.show()
