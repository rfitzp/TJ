# dT.py

# Plots delta T_e component of unreconnected eigenfunction associated with given rational surface in R, Z plane.
# User prompted for rational surface number and toroidal gridpoint.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

ncont = 160

fn   = 'Equilibrium1.nc'
ds   = nc.Dataset(fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']
t    = ds['theta']

epsa = para[0]
scale = 1.3*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'TJ1.nc'
ds1   = nc.Dataset(fn1)
rres  = ds1['r_res']
psi_r = ds1['delta_Te']
pn    = psi_r.shape[3]

fn2   = 'TJ2.nc'
ds2   = nc.Dataset(fn2)
rres1  = ds2['r_res']
psi_r1 = ds2['delta_Te']

fig = plt.figure (figsize = (12., 5.8))
fig.canvas.manager.set_window_title (r'TJ Code: delta T_e(R, Z)')
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17) 

nres = len(rres)

k = 0
n1 = 7

Pr0 = np.asarray(psi_r[k,:,:,n1])
Pr1 = np.asarray(psi_r1[k,:,:,n1])

plt.subplot (1, 2, 1)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres1, colors = 'black', linewidths = 0.5)    

cp = plt.contourf (RR, ZZ, Pr1, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "20")
plt.ylabel (r'$Z/R_0$', fontsize = "20")

plt.subplot (1, 2, 2)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

cp = plt.contourf (RR, ZZ, Pr0, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "20")
plt.ylabel (r'$Z/R_0$', fontsize = "20")

plt.tight_layout()

#plt.show()    
plt.savefig("dTea.png", dpi=300)
