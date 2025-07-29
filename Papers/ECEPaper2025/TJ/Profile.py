# Te.py

# Plots electron temperature and number density as well as gradients versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium1.nc'
ds = nc.Dataset(fn)
r  = ds['r']
q  = ds['q']
p2 = ds['p_2']
ne = ds['ne']
Te = ds['Te']

fn1  = 'TJ1.nc'
ds1  = nc.Dataset (fn1)
rres = np.asarray(ds1['r_res'])

fn2   = 'TJ2.nc'
ds2   = nc.Dataset (fn2)
rres2 = np.asarray(ds2['r_res'])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Electron Number Density and Temperature Profiles')
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20)
font = 20

plt.subplot (2, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (1.5,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (2.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (rres[0],  color = 'red', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (rres2[0], color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = font)
plt.ylabel (r'$q$',       fontsize = font)

plt.subplot (2, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, p2, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$q_0$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (rres[0],  color = 'red', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (rres2[0], color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = font)
plt.ylabel (r"$p_2$",     fontsize = font)

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, ne, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (rres[0],  color = 'red', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (rres2[0], color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',     fontsize = font)
plt.ylabel (r'$n_e(m^{-3})$', fontsize = font)

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (r, Te, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$q_0$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.axvline (rres[0],  color = 'red', linewidth = 1.5, linestyle = 'dashed')
plt.axvline (rres2[0], color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',  fontsize = font)
plt.ylabel (r"$T_e(eV)$",  fontsize = font)

plt.tight_layout ()

plt.savefig("Fig2.pdf")

#plt.show ()    
