# Lambda.py

# Plots eigenvalues of inverse plasma energy matrix versus r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
lam   = ds['lambda']
rres  = ds['r_res']

print (r.shape[0], mpol.shape[0], lam.shape[0], lam.shape[1])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Eigenvalues of Inverse Plasma Energy Matrix')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

plt.subplot(1, 1, 1)

plt.xlim(0., 1.)

plt.plot (r, lam[0,:], color = 'blue',    linewidth = 1, linestyle = 'solid', label = '$\lambda_0$')
plt.plot (r, lam[1,:], color = 'red',     linewidth = 1, linestyle = 'solid', label = '$\lambda_1$')
plt.plot (r, lam[2,:], color = 'green',   linewidth = 1, linestyle = 'solid', label = '$\lambda_2$')
plt.plot (r, lam[3,:], color = 'yellow',  linewidth = 1, linestyle = 'solid', label = '$\lambda_3$')
plt.plot (r, lam[4,:], color = 'magenta', linewidth = 1, linestyle = 'solid', label = '$\lambda_4$')
plt.plot (r, lam[5,:], color = 'cyan',    linewidth = 1, linestyle = 'solid', label = '$\lambda_5$')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$\lambda$", fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
