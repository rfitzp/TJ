# dT.py

# Plots delta T_e versus x and zeta
import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('prism')

ncont = 360

fn   = '../../Outputs/Island/Island.nc'
ds   = nc.Dataset(fn)
x    = np.asarray(ds['X'])
z    = np.asarray(ds['zeta']) /math.pi
t    = np.asarray(ds['delta_T'])

xx = np.sin(math.pi*z/2.)/2.

fn1   = '../../Outputs/TJ/TJ.nc'
ds1   = nc.Dataset(fn1)
rres  = ds1['r_res']
psi_r = ds1['dTe_cos']
psi_i = ds1['dTe_sin']

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: delta T_e(x, zeta)')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

plt.subplot (1, 1, 1)
plt.ylim (0., x[-1])
plt.xlim (0., 2.)

plt.plot (z, xx, color = 'black', linewidth = 3, linestyle = 'solid')

plt.contourf (x, z, t, ncont, cmap = ReBu)

plt.ylabel (r'$x/W$',       fontsize = "12")
plt.xlabel (r'$\zeta/\pi$', fontsize = "12")

plt.tight_layout()

plt.show()    
#plt.savefig("dTe.png")
