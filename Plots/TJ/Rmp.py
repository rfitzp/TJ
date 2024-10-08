# Rmp.py

# Plots poloidal harmonics of psi and Z components of ideal plasma response to rmp versus r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Psi_rmp_r']
Psi_i = ds['Psi_rmp_i']
Z_r   = ds['Z_rmp_r']
Z_i   = ds['Z_rmp_i']
rres  = ds['r_res']
mres  = ds['m_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Ideal Response to RMP')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = Psi_r.shape[0]

plt.subplot(2, 2, 1)

plt.xlim(0., 1.)

for jp in range (J):
    if (mpol[jp]) % 7 == 0:
        plt.plot (r, Psi_r[jp,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 1:
        plt.plot (r, Psi_r[jp,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 2:
        plt.plot (r, Psi_r[jp,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 3:
        plt.plot (r, Psi_r[jp,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 4:
        plt.plot (r, Psi_r[jp,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 5:
        plt.plot (r, Psi_r[jp,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 6:
        plt.plot (r, Psi_r[jp,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\psi_m$)", fontsize = "15")

plt.subplot(2, 2, 2)

plt.xlim(0., 1.)

for jp in range (J):
    if (mpol[jp]) % 7 == 0:
        plt.plot (r, Psi_i[jp,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 1:
        plt.plot (r, Psi_i[jp,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 2:
        plt.plot (r, Psi_i[jp,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 3:
        plt.plot (r, Psi_i[jp,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 4:
        plt.plot (r, Psi_i[jp,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 5:
        plt.plot (r, Psi_i[jp,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 6:
        plt.plot (r, Psi_i[jp,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Im($\psi_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]) % 7 == 0:
        plt.plot (r, Z_r[jp,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 1:
        plt.plot (r, Z_r[jp,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 2:
        plt.plot (r, Z_r[jp,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 3:
        plt.plot (r, Z_r[jp,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 4:
        plt.plot (r, Z_r[jp,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 5:
        plt.plot (r, Z_r[jp,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 6:
        plt.plot (r, Z_r[jp,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($Z_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim(0., 1.)

for jp in range (J):
    if (mpol[jp]) % 7 == 0:
        plt.plot (r, Z_i[jp,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 1:
        plt.plot (r, Z_i[jp,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 2:
        plt.plot (r, Z_i[jp,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 3:
        plt.plot (r, Z_i[jp,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 4:
        plt.plot( r, Z_i[jp,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 5:
        plt.plot (r, Z_i[jp,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]) % 7 == 6:
        plt.plot (r, Z_i[jp,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($Z_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
