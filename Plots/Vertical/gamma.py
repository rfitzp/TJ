# gamma.py

# Plots normalized resistive wall growth rate versus wall thickness parameter

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.optimize import brentq


fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
o     = ds['OutputParameters']

gammaw = o[1]

def f(x, d):

    if (d == 0):
        return x - gammaw
    else:
        return math.sqrt(x/d) * math.tanh (math.sqrt(x*d)) - gammaw

dd = np.linspace (0., 1., 100)

gg = []

for d in dd:
    g = brentq(lambda x: f(x, d), 0, 2*gammaw)
    gg.append (g)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Resistive Wall Mode Growth-Rate')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

plt.subplot(1, 1, 1)

plt.xlim (0., 1.)

plt.plot    (dd, gg, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\delta_w$',       fontsize = "15")
plt.ylabel (r"$\hat{\gamma}_w$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
    
    
