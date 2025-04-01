# bpert.py

# Plots perturbed magnetic field components of unreconnected eigenfunction associated
# with given rational surface versus theta at given radial gridpoint
# User prompted for rational surface number and gridpoint

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Equilibrium/Equilibrium.nc'
ds  = nc.Dataset(fn)

fn1 = '../../Outputs/TJ/TJ.nc'
ds1 = nc.Dataset(fn1)
t   = np.asarray(ds1['theta'])/math.pi
bRc = ds1['b_R_cos']
bRs = ds1['b_R_sin']
bZc = ds1['b_Z_cos']
bZs = ds1['b_Z_sin']
bPc = ds1['b_phi_cos']
bPs = ds1['b_phi_sin']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Components of perturbed magnetic field at plasma boundary')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

m = input ("rational surface number ? ")
k = int(m) - 1

j = bRc.shape[1]
g = input ("grid point (0, %4d) ? " % (j-1))

plt.subplot (3, 1, 1)

plt.xlim (0., 2.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (t, bRc[k,g,:], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (t, bRs[k,g,:], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta$',  fontsize = "15")
plt.ylabel (r'$b_R(T)$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim (0., 2.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (t, bZc[k,g,:], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (t, bZs[k,g,:], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta$',  fontsize = "15")
plt.ylabel (r'$b_Z(T)$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (0., 2.)

plt.plot    (t, bPc[k,g,:], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (t, bPs[k,g,:], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,             color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\theta$',     fontsize = "15")
plt.ylabel (r'$b_\phi(T)$',  fontsize = "15")
plt.legend (fontsize = "15")
             
plt.tight_layout ()

plt.show ()    
