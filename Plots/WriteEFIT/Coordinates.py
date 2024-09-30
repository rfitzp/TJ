# Coordinates.py

# Script plots EFIT r, w coordinate system in R, Z plane

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = '../../Outputs/Equilibrium/EFIT.nc'
ds = nc.Dataset(fn)

para = ds['RealParameters']
rr   = ds['r']
w    = ds['w']
c    = ds['cosw']
s    = ds['sinw']
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

fig = plt.figure (figsize = (2.0 * 3.5 /aspect ** 0.5, 2.0 * 3.5 * aspect ** 0.5))
fig.canvas.manager.set_window_title("WriteEFIT: r, w Coordinate System")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

XX, YY = np.meshgrid (r, z, indexing = 'ij')
RR     = np.asarray (rr)
WW     = np.asarray (w)
CW     = np.asarray (c)
SW     = np.asarray (s)

plt.subplot (2, 2, 1)

plt.contour (XX, YY, RR, levels = 40, colors = 'black', linewidths = 1, linestyles = 'solid')

plt.plot (rax, zax, 'ro', markersize = 2)

plt.plot (rb, zb, color = 'blue',  linewidth = 2)
plt.plot (rl, zl, color = 'black', linewidth = 4)
    
plt.xlabel ('$R$(m)',              fontsize = '15')
plt.ylabel ("$Z$(m)", labelpad=-9, fontsize = '15')
plt.title  ("$r$",                 fontsize = '14')

plt.subplot (2, 2, 2)

plt.contour (XX, YY, WW, levels = 80, colors = 'red', linewidths = 0.5, linestyles = 'solid')

plt.plot (rax, zax, 'ro', markersize = 2)

plt.plot (rb, zb, color = 'blue',  linewidth = 2)
plt.plot (rl, zl, color = 'black', linewidth = 4)
    
plt.xlabel ('$R$(m)',              fontsize = '15')
plt.ylabel ("$Z$(m)", labelpad=-9, fontsize = '15')
plt.title  ("$\omega$",            fontsize = '14')

plt.subplot (2, 2, 3)

plt.contour (XX, YY, CW, levels = 40, colors = 'black', linewidths = 1, linestyles = 'solid')

plt.plot (rax, zax, 'ro', markersize = 2)

plt.plot (rb, zb, color = 'blue',  linewidth = 2)
plt.plot (rl, zl, color = 'black', linewidth = 4)
    
plt.xlabel ('$R$(m)',              fontsize = '15')
plt.ylabel ("$Z$(m)", labelpad=-9, fontsize = '15')
plt.title  ("$\cos(\omega)$",      fontsize = '14')

plt.subplot (2, 2, 4)

plt.contour (XX, YY, SW, levels = 40, colors = 'black', linewidths = 1, linestyles = 'solid')

plt.plot (rax, zax, 'ro', markersize = 2)

plt.plot (rb, zb, color = 'blue',  linewidth = 2)
plt.plot (rl, zl, color = 'black', linewidth = 4)
    
plt.xlabel ('$R$(m)',              fontsize = '15')
plt.ylabel ("$Z$(m)", labelpad=-9, fontsize = '15')
plt.title  ("$\sin(\omega)$",      fontsize = '14')

plt.tight_layout(pad=0.5)

plt.show()

