# Ideale.py

# Plots poloidal harmonics of y and Z components of ideal eigenfunction versus r.
# User prompted for eigenfunction number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/Vertical/Vertical.nc'
ds    = nc.Dataset(fn)
r     = ds['r_grid']
Psi_r = ds['y_e_r']
Psi_i = ds['y_e_i']
Z_r   = ds['Z_e_r']
Z_i   = ds['Z_e_i']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Ideal Eigenfunction')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

J = Psi_r.shape[0]
print ("solution number = (%d .. %d)" % (0, J-1))
j = input ("solution number ? ")

plt.subplot(2, 2, 1)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp) % 7 == 0:
        plt.plot (r, Psi_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 1:
        plt.plot (r, Psi_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 2:
        plt.plot (r, Psi_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 3:
        plt.plot (r, Psi_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 4:
        plt.plot (r, Psi_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 5:
        plt.plot (r, Psi_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 6:
        plt.plot (r, Psi_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($y_m$)", fontsize = "15")

plt.subplot(2, 2, 2)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp) % 7 == 0:
        plt.plot (r, Psi_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 1:
        plt.plot (r, Psi_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 2:
        plt.plot (r, Psi_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 3:
        plt.plot (r, Psi_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 4:
        plt.plot (r, Psi_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 5:
        plt.plot (r, Psi_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 6:
        plt.plot (r, Psi_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($y_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

for jp in range (J):
    if (jp) % 7 == 0:
        plt.plot (r, Z_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 1:
        plt.plot (r, Z_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 2:
        plt.plot (r, Z_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 3:
        plt.plot (r, Z_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 4:
        plt.plot (r, Z_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 5:
        plt.plot (r, Z_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 6:
        plt.plot (r, Z_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($Z_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp) % 7 == 0:
        plt.plot (r, Z_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 1:
        plt.plot (r, Z_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 2:
        plt.plot (r, Z_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 3:
        plt.plot (r, Z_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 4:
        plt.plot( r, Z_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 5:
        plt.plot (r, Z_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (jp) % 7 == 6:
        plt.plot (r, Z_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($Z_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
