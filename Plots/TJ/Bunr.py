# Bunrc.py

# Plots poloidal harmonics of z and chi components of unreconnected solution vector associated with given rational surface versus r.
# User prompted for rational surface number.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['z_unrc_r']
Psi_i = ds['z_unrc_i']
Z_r   = ds['chi_unrc_r']
Z_i   = ds['chi_unrc_i']
rres  = ds['r_res']
mres  = ds['m_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Unreconnected Eigenfunctions')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize  =15)

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m = input ("rational surface number ? ")
j = int(m) - 1

J = Psi_r.shape[0]    

plt.subplot(2, 2, 1)

plt.xlim(0., 1.)

for jp in range (J):
    if (mpol[jp]-mres[j]) % 7 == 0:
        plt.plot (r, Psi_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 1:
        plt.plot (r, Psi_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 2:
        plt.plot (r, Psi_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 3:
        plt.plot (r, Psi_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 4:
        plt.plot (r, Psi_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 5:
        plt.plot (r, Psi_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 6:
        plt.plot (r, Psi_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15") 
plt.ylabel (r"Re($z_m$)", fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-mres[j]) % 7 == 0:
        plt.plot (r, Psi_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 1:
        plt.plot (r, Psi_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 2:
        plt.plot (r, Psi_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 3:
        plt.plot (r, Psi_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 4:
        plt.plot (r, Psi_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 5:
        plt.plot (r, Psi_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 6:
        plt.plot (r, Psi_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($z_m$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-mres[j]) % 7 == 0:
        plt.plot (r, Z_r[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 1:
        plt.plot (r, Z_r[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 2:
        plt.plot (r, Z_r[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 3:
        plt.plot (r, Z_r[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 4:
        plt.plot (r, Z_r[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 5:
        plt.plot (r, Z_r[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 6:
        plt.plot (r, Z_r[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\chi_m$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

for jp in range (J):
    if (mpol[jp]-mres[j]) % 7 == 0:
        plt.plot (r, Z_i[jp,j,:], color = 'black',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 1:
        plt.plot (r, Z_i[jp,j,:], color = 'red',     linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 2:
        plt.plot (r, Z_i[jp,j,:], color = 'green',   linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 3:
        plt.plot (r, Z_i[jp,j,:], color = 'blue',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 4:
        plt.plot (r, Z_i[jp,j,:], color = 'yellow',  linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 5:
        plt.plot (r, Z_i[jp,j,:], color = 'cyan',    linewidth = 1, linestyle = 'solid')
    elif (mpol[jp]-mres[j]) % 7 == 6:
        plt.plot (r, Z_i[jp,j,:], color = 'magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Im($\chi_m$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
