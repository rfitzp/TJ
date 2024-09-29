# Equilibrium2.py

# Script plots Stage2 plasma equilibrium gradients

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = '../../Outputs/Flux/Stage1.nc'
ds = nc.Dataset(fn)

para = ds['Parameters']
rl   = ds['RLIM']
zl   = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi  = para[5]
rax  = para[6]
zax  = para[7]

aspect = (zhi - zlow) / (rrgt - rlft)

fn1  = '../../Outputs/Flux/Stage2.nc'
ds1  = nc.Dataset(fn1)
psir = ds1['PSI_R']
psiz = ds1['PSI_Z']
r    = ds1['R']
z    = ds1['Z']
rb   = ds1['RBPTS']
zb   = ds1['ZBPTS']

fig = plt.figure (figsize=(2.0 * 6.0 /aspect ** 0.5, 6.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title ("FLUX: Equilibrium Flux-Surface Gradients")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

XX, YY = np.meshgrid (r, z, indexing = 'ij')
ZR     = np.asarray (psir)
zrmin  = np.min (ZR)
zrmax  = np.max (ZR)
levr   = np.linspace (zrmin, zrmax, 80)
ZZ     = np.asarray (psiz)
zzmin  = np.min (ZZ)
zzmax  = np.max (ZZ)
levz   = np.linspace (zzmin, zzmax, 80)
zlevel = np.array ([0.0])

plt.subplot (1, 2, 1)

plt.contour (XX, YY, ZR, levr, linewidths = 1.)
plt.contour (XX, YY, ZR, zlevel, colors = 'red', linewidths = 1.)
plt.plot (rax, zax,       'kx',    markersize = 4)
plt.plot (rb, zb, color = 'blue',  linewidth  = 1)
plt.plot (rl, zl, color = 'black', linewidth  = 4)

plt.xlabel ('$R/R_0$',  fontsize = '14')
plt.ylabel ("$Z/R_0$",  fontsize = '14')
plt.title ("$\\Psi_R$", fontsize = '14')

plt.subplot (1, 2, 2)

plt.contour (XX, YY, ZZ, levr, linewidths = 1.)
plt.contour (XX, YY, ZZ, zlevel, colors = 'red', linewidths = 1.)
plt.plot (rax, zax,       'kx',    markersize = 4)
plt.plot (rb, zb, color = 'blue',  linewidth  = 1)
plt.plot (rl, zl, color = 'black', linewidth  = 4)

plt.xlabel ('$R/R_0$',  fontsize = '14')
plt.ylabel ("$Z/R_0$",  fontsize = '14')
plt.title ("$\\Psi_Z$", fontsize = '14')

plt.tight_layout(pad=0.5)

plt.show()
