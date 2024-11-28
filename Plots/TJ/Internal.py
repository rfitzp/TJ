# Internal.py

# Plots poloidal harmonics of psi and Chi components of m-dominant internal ideal solution launched from magnetic axis versus r.
# User prompted for m.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Xi_i_r']
Psi_i = ds['Xi_i_i']
Chi_r = ds['Chi_i_r']
Chi_i = ds['Chi_i_i']
rres  = ds['r_res']
mres  = ds['m_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: mth Poloidal Harmonic Dominant Internal Ideal Solution Vector')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

print ("solution number = (%d .. %d)" % (mpol[0], mpol[-1]))
m   = input ("solution number ? ")
j   = int(m)  - mpol[0]

J = Psi_r.shape[0]

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-j) % 7 == 0:
        plt.plot (r, Psi_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 1:
        plt.plot (r, Psi_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 2:
        plt.plot (r, Psi_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 3:
        plt.plot (r, Psi_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 4:
        plt.plot (r, Psi_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 5:
        plt.plot (r, Psi_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 6:
        plt.plot (r, Psi_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',   fontsize = "15")
plt.ylabel (r"Re($\Xi_m$)", fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-j) % 7 == 0:
        plt.plot (r, Psi_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 1:
        plt.plot (r, Psi_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 2:
        plt.plot (r, Psi_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 3:
        plt.plot (r, Psi_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 4:
        plt.plot (r, Psi_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 5:
        plt.plot (r, Psi_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 6:
        plt.plot (r, Psi_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',            fontsize = "15")
plt.ylabel (r"Im($\Xi_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-j) % 7 == 0:
        plt.plot (r, Chi_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 1:
        plt.plot (r, Chi_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 2:
        plt.plot (r, Chi_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 3:
        plt.plot (r, Chi_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 4:
        plt.plot (r, Chi_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 5:
        plt.plot (r, Chi_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 6:
        plt.plot (r, Chi_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\chi_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-j) % 7 == 0:
        plt.plot (r, Chi_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 1:
        plt.plot (r, Chi_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 2:
        plt.plot (r, Chi_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 3:
        plt.plot (r, Chi_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 4:
        plt.plot( r, Chi_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 5:
        plt.plot (r, Chi_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-j) % 7 == 6:
        plt.plot (r, Chi_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Im($\chi_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
