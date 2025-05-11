# Te.py

# Plots electron temperature and number density as well as gradients versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = 'Equilibrium1.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
ne  = ds['ne']
nep = ds['nep']
Te  = ds['Te']
Tep = ds['Tep']
q   = ds['q']
p2  = ds['p_2']

fn1   = 'TJ1.nc'
ds1   = nc.Dataset (fn1)
rre1  = np.asarray(ds1['r_res'])

fn2   = 'TJ2.nc'
ds2   = nc.Dataset (fn2)
rre2  = np.asarray(ds2['r_res'])

rres = [rre2[0], rre1[0]]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Electron Number Density and Temperature Profiles')
plt.rc ('xtick', labelsize = 13) 
plt.rc ('ytick', labelsize = 13) 

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

#plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.5,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rr in rres:
    plt.axvline (rr, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',     fontsize = "15")
plt.ylabel (r'$q$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

#plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, p2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rr in rres:
    plt.axvline (rr, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$p_2$', fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, ne, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rr in rres:
    plt.axvline (rr, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',     fontsize = "15")
plt.ylabel (r'$n_e(m^{-3})$', fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, Te, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rr in rres:
    plt.axvline (rr, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$T_e(eV)$', fontsize = "15")

plt.tight_layout ()

#plt.show ()    
plt.savefig("Profile.pdf")
