# Profile2.py

# Script plots Stage2 equilibrium profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = '../../Outputs/Flux/Stage2.nc'
ds = nc.Dataset(fn)

r   = ds['r']
psi = ds['PsiN']
g   = ds['g']
P   = ds['P']
q   = ds['q']
f   = ds['f']

fig = plt.figure (figsize = (8.0, 6.0))
fig.canvas.manager.set_window_title ("FLUX: Stage2 Equilibrium Profiles")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

plt.subplot(2, 2, 1)

plt.xlim (0.0, 1.0)

plt.plot (psi, g, color = 'blue', linewidth = 1.5)

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ('$g$',       fontsize = 12)

plt.subplot (2, 2, 2)

plt.xlim (0.0, 1.0)

plt.plot (psi, P, color = 'blue', linewidth = 1.5)

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ('$P$',       fontsize = 12)

plt.subplot (2, 2, 3)

plt.xlim (0.0, 1.0)

plt.plot (psi, q, color = 'blue', linewidth = 1.5)

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$q$",       fontsize = 12)

plt.subplot (2, 2, 4)

plt.xlim (0.0, 1.0)

plt.plot (psi, f, color = 'blue', linewidth = 1.5)

plt.xlabel ('$\\Psi_N$', fontsize = 12)
plt.ylabel ("$f$",       fontsize = 12)

plt.tight_layout(pad=0.5)

plt.show()
