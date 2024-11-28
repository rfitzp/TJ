# Wall.py

# Plots data relating to wall

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Equilibrium/Equilibrium.nc'
ds = nc.Dataset(fn)
R  = ds['Rwall']
Z  = ds['Zwall']
w  = ds['wwall']
Rb = ds['Rbound']
Zb = ds['Zbound']
wb = ds['wbound']

ww  = np.asarray(w) /math.pi
RR  = np.asarray(R)  - 1.
wwb = np.asarray(wb) /math.pi
RRb = np.asarray(Rb) - 1. 

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Wall Data')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0., 2.)

plt.plot    (ww, RR,   color = 'blue',  linewidth = 2,    linestyle = 'solid',  label = r'$R_w-1$')
plt.plot    (ww, Z,    color = 'green', linewidth = 2,    linestyle = 'solid',  label = r'$Z_w$')
plt.plot    (wwb, RRb, color = 'blue',  linewidth = 2,    linestyle = 'dashed', label = r'$R_b-1$')
plt.plot    (wwb, Zb,  color = 'green', linewidth = 2,    linestyle = 'dashed', label = r'$Z_b$')
plt.axhline (0.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.,     color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (0.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')
plt.axvline (1.5,    color = 'black', linewidth = 1.5,  linestyle = 'dotted')

plt.xlabel (r'$\omega/\pi$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
