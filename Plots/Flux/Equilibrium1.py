# Equilibrium1.py

# Script plots Stage1 plasma equilibrium magnetic flux-surface in R, Z plane

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = '../../Outputs/Flux/Stage1.nc'
ds = nc.Dataset (fn)

para = ds['Parameters']
psi  = ds['PSI']
r    = ds['R']
z    = ds['Z']
rb   = ds['RBOUND']
zb   = ds['ZBOUND']
rl   = ds['RLIM']
zl   = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi  = para[5]
rax  = para[6]
zax  = para[7]
pax  = para[8]
pab  = para[9]

aspect = (zhi - zlow) / (rrgt - rlft)

fig = plt.figure (figsize = (7.0 /aspect ** 0.5, 7.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title ("FLUX: Equilibrium Magnetic Flux-Surfaces")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

XX, YY = np.meshgrid (r, z, indexing = 'ij')
ZZ     = np.asarray (psi)

plt.contour (XX, YY, ZZ, 40, colors = 'black', linewidths = 0.5, linestyles = 'solid')

plt.plot (rax, zax, 'ro',            markersize = 2)
plt.plot (rb, zb,   color = 'blue',  linewidth  = 2)
plt.plot (rl, zl,   color = 'black', linewidth  = 4)
    
plt.xlabel ('$R/R_0$', fontsize = '15')
plt.ylabel ("$Z/R_0$", labelpad = -9, fontsize = '15')

plt.show()

