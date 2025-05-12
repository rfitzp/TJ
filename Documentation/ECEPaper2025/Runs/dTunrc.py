# dTunr.py

# Plots poloidal harmonics of delta T_e and delta n_e  components of unreconnected solution vector associated with given rational surface versus r.
# User prompted for rational surface number.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ1.nc'
ds   = nc.Dataset(fn)
mpol = ds['mpol']
r    = ds['r_grid']
Te_r = ds['delta_Te_unrc_r']
Te_i = ds['delta_Te_unrc_i']
ne_r = ds['delta_ne_unrc_r']
ne_i = ds['delta_ne_unrc_i']
rres = ds['r_res']
mres = ds['m_res']

fn1  = 'TJ2.nc'
ds1   = nc.Dataset(fn1)
mpol1 = ds1['mpol']
r1    = ds1['r_grid']
Te_r1 = ds1['delta_Te_unrc_r']
Te_i1 = ds1['delta_Te_unrc_i']
ne_r1 = ds1['delta_ne_unrc_r']
ne_i1 = ds1['delta_ne_unrc_i']
rres1 = ds1['r_res']
mres1 = ds1['m_res']

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Unreconnected Eigenfunctions')
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

j = 0

plt.subplot(2, 1, 1)

plt.xlim(0., 1.)

plt.plot (r, Te_r1[9,j,:],  color = 'cyan',    linewidth = 2, linestyle = 'solid', label = "$m=-1$")
plt.plot (r, Te_r1[11,j,:], color = 'magenta', linewidth = 2, linestyle = 'solid', label = "$m=1$")
plt.plot (r, Te_r1[12,j,:], color = 'red',     linewidth = 2, linestyle = 'solid', label = "$m=2$")
plt.plot (r, Te_r1[13,j,:], color = 'black',   linewidth = 2, linestyle = 'solid', label = "$m=3$")
plt.plot (r, Te_r1[14,j,:], color = 'green',   linewidth = 2, linestyle = 'solid', label = "$m=4$")
plt.plot (r, Te_r1[15,j,:], color = 'blue',    linewidth = 2, linestyle = 'solid', label = "$m=5$")
plt.plot (r, Te_r1[16,j,:], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = "$m=6$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres1:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.ylabel (r'$T_e(eV)$', fontsize = fontsize) 
plt.legend (fontsize = 12)

plt.subplot(2, 1, 2)

plt.xlim(0., 1.)

plt.plot (r, Te_r[8,j,:],  color = 'cyan',    linewidth = 2, linestyle = 'solid', label = "$m=-2$")
plt.plot (r, Te_r[9,j,:],  color = 'magenta', linewidth = 2, linestyle = 'solid', label = "$m=-1$")
plt.plot (r, Te_r[11,j,:], color = 'red',     linewidth = 2, linestyle = 'solid', label = "$m=1$")
plt.plot (r, Te_r[12,j,:], color = 'black',   linewidth = 2, linestyle = 'solid', label = "$m=2$")
plt.plot (r, Te_r[13,j,:], color = 'green',   linewidth = 2, linestyle = 'solid', label = "$m=3$")
plt.plot (r, Te_r[14,j,:], color = 'blue',    linewidth = 2, linestyle = 'solid', label = "$m=4$")
plt.plot (r, Te_r[15,j,:], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = "$m=5$")

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.ylabel (r'$T_e(eV)$', fontsize = fontsize) 
plt.legend (fontsize = 12)

plt.tight_layout ()

#plt.show ()    
plt.savefig ("dTunrc.pdf")
