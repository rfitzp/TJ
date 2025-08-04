# alpha.py

# Plots quantities as function of alpha

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/StartUp/StartUp.nc'
ds = nc.Dataset (fn)
a  = ds['alpha']
b  = ds['lambda']
q  = ds['q_a']
l  = 2*np.asarray(ds['l_i'])
T  = ds['T_ramp']
t  = ds['t_ramp']
E  = ds['E_ramp']
P  = ds['P_ramp']

font = 17
fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'StartUp Code: alpha scan')
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 1, 1)

plt.xlim (0., a[-1])
 
plt.plot    (a, q, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = '$q_a$')
plt.plot    (a, b, color = 'red',   linewidth = 2,  linestyle = 'solid', label = '$\lambda$')

plt.xlabel (r'$\alpha$', fontsize = font)
plt.legend (fontsize = 15, loc = "upper left", ncol = 2)

plt.subplot (3, 1, 2)

plt.xlim (0., a[-1])
 
plt.plot    (a, l, color = 'blue',    linewidth = 2, linestyle = 'solid',  label = '$2\,l_i$')
plt.plot    (a, T, color = 'red',     linewidth = 2, linestyle = 'solid',  label = '$T_{ramp}$')
plt.plot    (a, E, color = 'green',   linewidth = 2, linestyle = 'solid',  label = '${\cal E}_{ramp}$')

plt.xlabel (r'$\alpha$', fontsize = font)
plt.legend (fontsize = 15, loc = "upper left", ncol = 3)

plt.subplot (3, 1, 3)

plt.xlim (0., a[-1])
 
plt.plot    (a, t, color = 'blue',  linewidth = 2, linestyle = 'solid')

plt.xlabel (r'$\alpha$',            fontsize = font)
plt.ylabel (r"$\tau_{min}/\tau_R$", fontsize = font)

plt.tight_layout ()

plt.show () 
#plt.savefig("Figure2.pdf")
