# Profile1.py

# Script plots Stage1 plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn  = '../../Outputs/Flux/Stage1.nc'
ds  = nc.Dataset(fn)

psin = ds['PSI_N']
p    = ds['p']
pp   = ds['pp']
t    = ds['t']
tp   = ds['ttp']
q    = ds['q']

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title ("FLUX: Stage1 Equilibrium Profiles")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

plt.subplot (2, 3, 1)

plt.xlim (0.0, 1.0)

plt.plot (psin, p, color = 'blue', linewidth = 1.5)
plt.axhline (0.,   color = 'red',  linewidth = 1, linestyle = 'dotted')

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ('$P$',       fontsize = 12)

plt.subplot(2, 3, 2)

plt.xlim (0.0, 1.0)

plt.plot (psin, t, color = 'blue', linewidth = 1.5)
plt.axhline (1.,   color = 'red',  linewidth = 1, linestyle='dotted')

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$g$",       fontsize = 12)

plt.subplot (2, 3, 3)

plt.xlim (0.0, 1.0)

plt.plot (psin, q, color = 'blue', linewidth = 1.5)

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$q$",       fontsize = 12)

plt.subplot (2, 3, 4)

plt.xlim (0.0, 1.0)

plt.plot (psin, pp, color = 'blue', linewidth = 1.5)
plt.axhline (0.,    color = 'red',  linewidth = 1, linestyle='dotted')

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$P'$",      fontsize = 12)

plt.subplot (2, 3, 5)

plt.xlim (0.0, 1.0)

plt.plot (psin, tp, color = 'blue', linewidth = 1.5)
plt.axhline (0.,    color = 'red',  linewidth = 1, linestyle='dotted')

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$gg'$",     fontsize = 12)

plt.tight_layout(pad=0.5)

plt.show()
