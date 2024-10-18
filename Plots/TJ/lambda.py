# Lambda.py

# Plots eigenvalues of inverse plasma energy matrix versus r.
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
fig.canvas.manager.set_window_title (r'TJ Code: Eigenvalues of Plasma Energy Matrix')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

llam = np.asarray(lam)

plt.xlim (0., 1.)

plt.ylim (1.1*min(llam[0,:]), 1.1*llam[-1,-1])

plt.plot (r, lam[0, :],  color = 'blue', linewidth = 1, linestyle = 'solid', label = '$\lambda_{min}$')
plt.plot (r, lam[-1, :], color = 'red',  linewidth = 1, linestyle = 'solid', label = '$\lambda_{max}$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$\lambda$", fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)

for j in range(1,J-1):
    if j % 7 == 0:
        plt.plot (r, lam[j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif j % 7 == 1:
        plt.plot (r, lam[j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif j % 7 == 2:
        plt.plot (r, lam[j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif j % 7 == 3:
        plt.plot (r, lam[j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif j % 7 == 4:
        plt.plot( r, lam[j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif j % 7 == 5:
        plt.plot (r, lam[j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif j % 7 == 6:
        plt.plot (r, lam[j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$\lambda$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
