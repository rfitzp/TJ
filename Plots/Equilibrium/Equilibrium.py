# Equilibrium.py

# Plots components of aspect-ratio expanded equilibrium versus r

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Equilibrium/Equilibrium.nc'
ds = nc.Dataset(fn)
r  = ds['r']
g2 = ds['g_2']
p2 = ds['p_2']
f1 = ds['f_1']
f3 = ds['f_3']
q0 = ds['q_0']
q2 = ds['q_2']
ip = ds['I_p']
it = ds['I_t']
jp = ds['J_p']
jt = ds['J_t']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Equilibrium Quantities')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (4, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, p2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$p_2$',     fontsize = "15")

plt.subplot (4, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, g2, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$q_0$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$g_2$',     fontsize = "15")

plt.subplot (4, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r, f1, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$f_1$")
plt.plot    (r, f3, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$f_3$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.legend (fontsize = '15')

plt.subplot (4, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r, q0, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$q_0$")
plt.plot    (r, q2, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$q_2$")
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.legend (fontsize = '15')

plt.subplot (4, 2, 5)

plt.xlim (0., 1.)

plt.plot    (r, it, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r'$\hat{I}_t$',  fontsize = "15")

plt.subplot (4, 2, 6)

plt.xlim (0., 1.)

plt.plot    (r, ip, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',   fontsize = "15")
plt.ylabel (r'$\hat{I}_p$', fontsize = "15")

plt.subplot (4, 2, 7)

plt.xlim (0., 1.)

plt.plot    (r, jt, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"$\hat{I}_t'$", fontsize = "15")

plt.subplot (4, 2, 8)

plt.xlim (0., 1.)

plt.plot    (r, jp, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r"$\hat{I}_p'$", fontsize = "15")

plt.tight_layout ()
plt.show ()    
