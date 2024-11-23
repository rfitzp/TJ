# Wall.py

# Plots data relating to wall

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Equilibrium/Equilibrium.nc'
ds   = nc.Dataset(fn)
R    = ds['Rwall']
Z    = ds['Zwall']
w    = ds['wwall']

ww = np.asarray(w) /math.pi
RR = np.asarray(R) - 1.  

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Wall Data')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0., 2.)

plt.plot    (ww, RR, color = 'blue',  linewidth = 2,    linestyle = 'solid', label = r'$R-1$')
plt.plot    (ww, Z,  color = 'green', linewidth = 2,    linestyle = 'solid', label = r'$Z$')
plt.axhline (0.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')

plt.xlabel (r'$\omega/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
