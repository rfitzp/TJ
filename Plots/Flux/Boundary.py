# Boundary.py

# Plots data relating to plasma/vacuum interface

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Flux/Stage2.nc'
ds   = nc.Dataset(fn)
t    = ds['tbound']
R    = ds['Rbound']
Z    = ds['Zbound']
Rt   = ds['dRdtheta']
Zt   = ds['dZdtheta']

fn   = '../../Outputs/Equilibrium/Equilibrium.nc'
ds   = nc.Dataset(fn)
t1    = ds['tbound']
Rt1   = ds['dRdtheta']
Zt1   = ds['dZdtheta']

tt  = np.asarray(t)  /math.pi
tt1 = np.asarray(t1) /math.pi
RR  = np.asarray(R) - 1.  

fig = plt.figure (figsize = (12.0, 5.0))
fig.canvas.manager.set_window_title (r'Flux: Plasma/Vacuum Interface Data')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 2, 1)

plt.xlim (0., 2.)

plt.plot    (tt, RR, color = 'blue',  linewidth = 2,    linestyle = 'solid', label = r'$R-1$')
plt.plot    (tt, Z,  color = 'green', linewidth = 2,    linestyle = 'solid', label = r'$Z$')
plt.axhline (0.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (1, 2, 2)

plt.xlim (0., 2.)

plt.plot    (tt, Rt,   color = 'blue',  linewidth = 1,   linestyle = 'solid', label = r"$\partial R/\partial\theta$")
plt.plot    (tt, Zt,   color = 'green', linewidth = 1,   linestyle = 'solid', label = r"$\partial Z/\partial\theta$")
plt.plot    (tt1, Rt1, color = 'blue',  linewidth = 1,   linestyle = 'dotted')
plt.plot    (tt1, Zt1, color = 'green', linewidth = 1,   linestyle = 'dotted')
plt.axvline (1.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.5,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
