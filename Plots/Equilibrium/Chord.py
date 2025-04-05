# Chord.py

# Plots equilibrium quantities along tilted central chord

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn1 = '../../Outputs/Equilibrium/Equilibrium.nc'
ds1 = nc.Dataset(fn1)
req  = ds1['r_eq']
weq  = ds1['omega_eq']
teq  = ds1['theta_eq']
Req  = ds1['R_eq']
Zeq  = ds1['Z_eq']
BReq = ds1['BR_eq']
neeq = ds1['ne_eq']
Teeq = ds1['Te_eq']

RReq = np.asarray(Req) - 1.
wweq = np.asarray(weq) /math.pi
tteq = np.asarray(teq) /math.pi

Nf = len(RReq)

x = np.linspace (0., 1., Nf);

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Equilibrium Quantities Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (0., 1.)

plt.plot    (x, req, color = 'blue',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,    color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',       fontsize = "15")
plt.ylabel (r'$\hat{r}$', fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (0., 1.)

plt.plot    (x, wweq, color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r'$\omega/\pi$')
plt.plot    (x, tteq, color = 'red',   linewidth = 2,   linestyle = 'solid',  label = r'$\theta/\pi$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 1.)

plt.plot    (x, RReq, color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r'$R-1$')
plt.plot    (x, Zeq,  color = 'red',   linewidth = 2,   linestyle = 'solid',  label = r'$Z$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 1.)

plt.plot    (x, BReq, color = 'blue',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',              fontsize = "15")
plt.ylabel (r'$B_\parallel(T)$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (x, neeq, color = 'blue',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',      fontsize = "15")
plt.ylabel (r'$n_e(m^{-3})$', fontsize = "15")

plt.subplot (3, 2, 6)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (x, Teeq, color = 'blue',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$',       fontsize = "15")
plt.ylabel (r'$T_e(eV)$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
