# Chord.py

# Plots ECE frequencies along tilted central chord

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn1 = '../../Outputs/Equilibrium/Equilibrium.nc'
ds1 = nc.Dataset(fn1)
we  = ds1['We_eq']
wu  = ds1['wU_eq']
wl  = ds1['wL_eq']
wuh = ds1['wUH_eq']

we2 = 2.*np.asarray(we)

Nf = len(we)

x = np.linspace (0., 1., Nf);

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: ECE Frequencies Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)

plt.plot    (x, we2, color = 'blue',   linewidth = 2,   linestyle = 'solid', label = r"$2\,\Omega_e$")
plt.plot    (x, we,  color = 'cyan',   linewidth = 2,   linestyle = 'solid', label = r"$\Omega_e$")
plt.plot    (x, wu,  color = 'green',  linewidth = 2,   linestyle = 'solid', label = r"$\omega_1$")
plt.plot    (x, wl,  color = 'red',    linewidth = 2,   linestyle = 'solid', label = r"$\omega_2$")
plt.plot    (x, wuh, color = 'black',  linewidth = 2,   linestyle = 'solid', label = r"$\omega_{UH}$")

plt.axhline (0.,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.,     color = 'black',  linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.5,    color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
