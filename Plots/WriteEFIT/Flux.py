# Flux.py

# Script plots EFIT magnetic flux-surfaces in R, Z plane

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = '../../Outputs/Equilibrium/EFIT.nc'
ds = nc.Dataset(fn)

para = ds['RealParameters']
psi  = ds['PSI']
r    = ds['RGRID']
z    = ds['ZGRID']
rb   = ds['RBOUND']
zb   = ds['ZBOUND']
rl   = ds['RLIMITER']
zl   = ds['ZLIMITER']
rlft = para[11]
rrgt = para[12]
zlow = para[13]
zhi  = para[14]
rax  = para[6]
zax  = para[7]
pax  = para[8]
pab  = para[9]

aspect = (zhi - zlow) / (rrgt - rlft)

fig = plt.figure (figsize = (7.0 /aspect ** 0.5, 7.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title("WriteEFIT: Magnetic Flux-Surfaces")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

XX, YY = np.meshgrid (r, z, indexing = 'ij')
ZZ     = np.asarray  (psi)

plt.contour (XX, YY, ZZ, levels = 40, colors = 'black', linewidths = 1, linestyles = 'solid')

plt.plot (rax, zax, 'ro', markersize = 2)

plt.plot (rb, zb, color = 'blue',  linewidth = 2)
plt.plot (rl, zl, color = 'black', linewidth = 4)
    
plt.xlabel ('$R$(m)',              fontsize = '15')
plt.ylabel ("$Z$(m)", labelpad=-9, fontsize = '15')

plt.show()

