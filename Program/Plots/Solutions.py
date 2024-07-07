# Solutions.py

# Plots psi and Z components of m-dominant solution vector versus r.
# User prompted for m.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['Psi_r']
Psi_i = ds['Psi_i']
Z_r   = ds['Z_r']
Z_i   = ds['Z_i']
Pn    = ds['Psi_norm']
Zn    = ds['Z_norm']
Tq    = ds['Torque_test']
rres  = ds['rres']

fig = plt.figure (figsize=(12.0, 8.0))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15)

nres = len(rres)
print ("solution number = (%d .. %d) or (%d .. %d)" % (mpol[0], mpol[-1], mpol[-1]+1, mpol[-1]+nres))
m   = input ("solution number ? ")
j   = int(m)  - mpol[0]

J = Psi_r.shape[0]    

plt.subplot(3, 2, 1)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp-j) % 7 == 0:
        plt.plot(r, Psi_r[jp,j,:], color='black', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 1:
        plt.plot(r, Psi_r[jp,j,:], color='red', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 2:
        plt.plot(r, Psi_r[jp,j,:], color='green', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 3:
        plt.plot(r, Psi_r[jp,j,:], color='blue', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 4:
        plt.plot(r, Psi_r[jp,j,:], color='yellow', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 5:
        plt.plot(r, Psi_r[jp,j,:], color='cyan', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 6:
        plt.plot(r, Psi_r[jp,j,:], color='magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"real($\psi$)", fontsize="15")

plt.subplot(3, 2, 2)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp-j) % 7 == 0:
        plt.plot(r, Psi_i[jp,j,:], color='black', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 1:
        plt.plot(r, Psi_i[jp,j,:], color='red', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 2:
        plt.plot(r, Psi_i[jp,j,:], color='green', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 3:
        plt.plot(r, Psi_i[jp,j,:], color='blue', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 4:
        plt.plot(r, Psi_i[jp,j,:], color='yellow', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 5:
        plt.plot(r, Psi_i[jp,j,:], color='cyan', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 6:
        plt.plot(r, Psi_i[jp,j,:], color='magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"imag($\psi$)", fontsize="15")

plt.subplot(3, 2, 3)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp-j) % 7 == 0:
        plt.plot(r, Z_r[jp,j,:], color='black', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 1:
        plt.plot(r, Z_r[jp,j,:], color='red', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 2:
        plt.plot(r, Z_r[jp,j,:], color='green', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 3:
        plt.plot(r, Z_r[jp,j,:], color='blue', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 4:
        plt.plot(r, Z_r[jp,j,:], color='yellow', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 5:
        plt.plot(r, Z_r[jp,j,:], color='cyan', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 6:
        plt.plot(r, Z_r[jp,j,:], color='magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"real($Z$)", fontsize="15")

plt.subplot(3, 2, 4)

plt.xlim(0., 1.)

for jp in range (J):
    if (jp-j) % 7 == 0:
        plt.plot(r, Z_i[jp,j,:], color='black', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 1:
        plt.plot(r, Z_i[jp,j,:], color='red', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 2:
        plt.plot(r, Z_i[jp,j,:], color='green', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 3:
        plt.plot(r, Z_i[jp,j,:], color='blue', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 4:
        plt.plot(r, Z_i[jp,j,:], color='yellow', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 5:
        plt.plot(r, Z_i[jp,j,:], color='cyan', linewidth = 1, linestyle = 'solid')
    elif (jp-j) % 7 == 6:
        plt.plot(r, Z_i[jp,j,:], color='magenta', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"imag($Z$)", fontsize="15")

plt.subplot(3, 2, 5)

plt.xlim(0., 1.)

plt.plot(r, Pn[j,:], color='blue', linewidth = 1, linestyle = 'solid')
plt.plot(r, Zn[j,:], color='red', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"norm($\psi$)", fontsize="15")

plt.subplot(3, 2, 6)

plt.xlim(0., 1.)

plt.plot(r, Tq[j,:], color='blue', linewidth = 1, linestyle = 'solid')

plt.axhline (0., color='black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color='black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel(r'$\hat{r}$', fontsize="15")
plt.ylabel(r"$T_\phi$", fontsize="15")
    
plt.tight_layout()

plt.show()    
