# alpha.py

# Plots quantities as function of alpha

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Startup.nc'
ds = nc.Dataset (fn)
a  = ds['alpha']
b  = ds['lambda']
q  = ds['q_a']
l  = ds['l_i']
T  = ds['T_ramp']
t  = ds['t_ramp']
E  = ds['E_ramp']

font = 20
fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'StartUp Code: alpha scan')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 1, 1)

plt.xlim (0., a[-1])
 
plt.plot    (a, q, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = '$q_a$')
plt.plot    (a, b, color = 'red',   linewidth = 2,  linestyle = 'solid', label = '$\lambda$')
plt.axhline (0.,   color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axhline (2.,   color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axhline (4.,   color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axhline (6.,   color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axhline (8.,   color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axhline (10.,  color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\alpha$', fontsize = font)
plt.legend (fontsize = 15)

plt.subplot (3, 1, 2)

plt.xlim (0., a[-1])
 
plt.plot    (a, l, color = 'blue',  linewidth = 2, linestyle = 'solid', label = '$l_i$')
plt.plot    (a, T, color = 'red',   linewidth = 2, linestyle = 'solid', label = '$T_{ramp}$')
plt.plot    (a, E, color = 'green', linewidth = 2, linestyle = 'solid', label = '$E_{ramp}$')
plt.axhline (0.,   color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$\alpha$', fontsize = font)
plt.legend (fontsize = 15)

plt.subplot (3, 1, 3)

plt.xlim (0., a[-1])
 
plt.plot    (a, t, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$\alpha$',  fontsize = font)
plt.ylabel ("$t_{ramp}$", fontsize = font)

plt.tight_layout ()

plt.show () 
