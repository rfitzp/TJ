# crit.py

# Plots lowest eigenvalue of inverse plasma energy matrix versus r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
r     = ds['r_grid']
lam   = ds['lambda']
rres  = ds['r_res']
mpol  = ds['mpol']

J = mpol.shape[0]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Lowest Eigenvalue of Plasma Energy Matrix')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

llam = np.asarray(lam)

plt.xlim (0., 1.)

plt.plot (r, lam[0, :],  color = 'blue', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',   fontsize = "15")
plt.ylabel (r"$\lambda_0$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
