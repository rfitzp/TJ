# dT.py

# Plots delta T_e versus x and zeta

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic') #plt.get_cmap ('prism')

ncont = 360

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset (fn)
x  = np.asarray (ds['X'])
z  = np.asarray (ds['zeta']) /math.pi
t  = np.asarray (ds['delta_T'])
p  = np.asarray (ds['InputParameters'])

delta = p[5];

xp = []
for y in z:
    val = (delta/8.**0.5) * math.cos(y*math.pi) + math.sqrt(1. - math.cos(y*math.pi - delta*delta*math.sin(y*math.pi))) /8.**0.5
    xp.append (val)

xm = []
for y in z:
    val = (delta/8.**0.5) * math.cos(y*math.pi) - math.sqrt(1. - math.cos(y*math.pi - delta*delta*math.sin(y*math.pi))) /8.**0.5
    xm.append (val)

x0 = []
for y in z:
    x0.append (0.)

fig = plt.figure (figsize = (12.0, 7.0))
fig.canvas.manager.set_window_title (r'Island Code: delta T_e(zeta, x)')
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17) 

plt.subplot (1, 1, 1)
plt.ylim (-1., 1.)
plt.xlim (z[0], z[-1])

plt.plot (z, xp, color = 'black', linewidth = 3, linestyle = 'solid')
plt.plot (z, xm, color = 'black', linewidth = 3, linestyle = 'solid')
plt.plot (z, x0, color = 'black', linewidth = 2, linestyle = 'dotted')

plt.plot ([1.], [-delta/8.**0.5], marker = 'o', markersize = 4, color = "black")

plt.contourf (z, x, t, ncont, cmap = ReBu)

levsp = np.arange (0.01, 1., 0.1)
levsm = - np.arange (0.02, 1., 0.1)
levs  = np.concatenate((levsm[::-1], levsp))

plt.contour (z, x, t, levels = levs, colors = 'black', linestyles = 'solid', linewidths = 0.5)  

plt.ylabel (r'$x/W$',       fontsize = "17")
plt.xlabel (r'$\zeta/\pi$', fontsize = "17")

plt.tight_layout ()

plt.show ()    
#plt.savefig("dTe.png")
