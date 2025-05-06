# Btoroidal.py

# Plots B_toroidal and |B| in R, Z plane.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

ncont = 160

fn   = '../../Outputs/Equilibrium/Equilibrium.nc'
ds   = nc.Dataset(fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']
t    = ds['theta']

epsa = para[0]
scale = 1.5*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = '../../Outputs/TJ/TJ.nc'
ds1   = nc.Dataset(fn1)
rres  = ds1['r_res']
btor  = np.asarray(ds1['B_toroidal'])
bmod  = np.asarray(ds1['B_modulus'])

fig = plt.figure (figsize = (12.5, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: B_toroidal (R, Z) and |B| (R, Z)')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

plt.subplot (1, 2, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

cp = plt.contourf (RR, ZZ, btor, ncont, cmap = ReBu)

plt.contour (RR, ZZ, bmod, rres, levels = 20, colors = 'white', linewidths = 0.5)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot (1, 2, 2)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)

cp = plt.contourf (RR, ZZ, bmod, ncont, cmap = ReBu)

plt.contour (RR, ZZ, bmod, rres, levels = 20, colors = 'white', linewidths = 0.5)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout()

plt.show()    
