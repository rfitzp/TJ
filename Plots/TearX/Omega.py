# Omega.py

# Plots diamagnetic and ExB frequencies versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = '../../Outputs/TearX/TearX.nc'
ds     = nc.Dataset(fn)
r      = ds['r']
omegae = ds['omegae']
omegai = ds['omegai']
omegaE = ds['omegaE']

oe = np.asarray(omegae) + np.asarray(omegaE)

rres = ds['rres']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: Frequencies')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)

plt.xlim (0., 1.)
 
plt.plot (r, omegae, color = 'red',   linewidth = 2,    linestyle = 'solid', label = r"$\omega_{\ast\,e}$")
plt.plot (r, omegai, color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r"$\omega_{\ast\,i}$")
plt.plot (r, omegaE, color = 'green', linewidth = 2,   linestyle = 'solid',  label = r"$\omega_E$")
plt.plot (r, oe,     color = 'cyan',  linewidth = 2,   linestyle = 'solid',  label = r"$\omega_{\perp\,e}$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
