# dTunr1.py

# Plots kth poloidal harmonic of delta T_e and delta n_e components of unreconnected solution vector associated with given rational surface versus r.
# User prompted for rational surface number and poloidal harmomic mode number.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = '../../Outputs/TJ/TJ.nc'
ds    = nc.Dataset(fn)
mpol  = ds['mpol']
r     = ds['r_grid']
Psi_r = ds['delta_Te_unrc_r']
Psi_i = ds['delta_Te_unrc_i']
Z_r   = ds['delta_ne_unrc_r']
Z_i   = ds['delta_ne_unrc_i']
rres  = ds['r_res']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Unreconnected Eigenfunctions')
plt.rc('xtick', labelsize = 15) 
plt.rc('ytick', labelsize = 15)

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m   = input ("rational surface number ? ")
j   = int(m) - 1
print ("poloidal mode number = (%d .. %d) (%d is 0/0)" % (mpol[0], mpol[-1]+1, mpol[-1]+1))
mp   = input ("m  ? ")
jp   = int(mp) - mpol[0]

J = Psi_r.shape[0]    

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot (r, Psi_r[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Re($\delta T_{e\,m}$)", fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim(0., 1.)

plt.plot(r, Psi_i[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"Im($\delta T_{e\,m}$)", fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot (r, Z_r[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Re($\delta n_{e\,m}$)", fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot (r, Z_i[jp,j,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"Im($\delta n_{e\,m}$)", fontsize = "15")

plt.tight_layout ()

plt.show ()    
