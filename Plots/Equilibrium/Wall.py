# Wall.py

# Plots data relating to wall

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Equilibrium/Equilibrium.nc'
ds   = nc.Dataset(fn)
t    = ds['twall']
R    = ds['Rwall']
Z    = ds['Zwall']
w    = ds['wwall']
R2   = ds['R2wall']
grr2 = ds['grr2wall']
Rt   = ds['dRdthetawall']
Zt   = ds['dZdthetawall']

tt = np.asarray(t) /math.pi
ww = np.asarray(w) /math.pi
RR = np.asarray(R) - 1.  

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Wall Data')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0., 2.)

plt.plot    (tt, RR, color = 'blue',  linewidth = 2,    linestyle = 'solid', label = r'$R-1$')
plt.plot    (tt, Z,  color = 'green', linewidth = 2,    linestyle = 'solid', label = r'$Z$')
plt.axhline (0.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 2.)
plt.ylim (0., 2.)

plt.plot    (ww, tt, color = 'blue',  linewidth = 2,    linestyle = 'solid')
plt.axhline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axhline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axhline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')

plt.plot ([0., 2.], [0., 2.], color = 'red', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\omega/\pi$', fontsize = "15")
plt.ylabel (r'$\theta/\pi$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 2.)

plt.plot    (tt, Rt, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$\partial R/\partial\theta$")
plt.plot    (tt, Zt, color = 'green', linewidth = 2,   linestyle = 'solid', label = r"$\partial Z/\partial\theta$")
plt.axvline (1.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 2.)

plt.plot    (tt, R2,   color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r'$R^2$')
plt.plot    (tt, grr2, color = 'green', linewidth = 2,   linestyle = 'solid', label = r'$|\nabla r|^2$')
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.5,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
