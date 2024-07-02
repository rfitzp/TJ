# RMP.py

# Plots magnetic scalar potential of resonant magnetic perturbation associated with rational surface in R, Z plane
# User prompted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap('seismic')

ncont = 160

fn   = 'Equilibrium.nc'
ds   = nc.Dataset(fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']
t    = ds['theta']

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1  = 'TJ.nc'
ds1  = nc.Dataset(fn1)
rres = ds1['rres']
Rv   = ds1['Rv']
Zv   = ds1['Zv']
Vxr  = ds1['Vx_r']
Vxi  = ds1['Vx_i']

rv   = np.asarray(Rv)
zv   = np.asarray(Zv)

[RRV, ZZV] = np.meshgrid(rv, zv)

fig = plt.figure(figsize=(12.0, 6.0))
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 

m = input ("rational surface number ? ")
k = int(m) - 1

Vr = Vxr[k,:,:]
Vi = Vxi[k,:,:]

plt.subplot(1, 2, 1)
plt.xlim(Rv[0], Rv[-1])
plt.ylim(Zv[0], Zv[-1])

plt.plot(R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour(RR, ZZ, rr, rres, colors='black', linewidths = 0.5)

plt.contourf(RRV, ZZV, Vr, ncont, cmap=ReBu)    

plt.plot([1.], [0.], marker='o', markersize=1, color="black")

plt.xlabel(r'$R/R_0$', fontsize="12")
plt.ylabel(r'$Z/R_0$',  fontsize="12")

plt.subplot(1, 2, 2)
plt.xlim(Rv[0], Rv[-1])
plt.ylim(Zv[0], Zv[-1])

plt.plot(R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour(RR, ZZ, rr, rres, colors='black', linewidths = 0.5)

plt.contourf(RRV, ZZV, Vi, ncont, cmap=ReBu)    

plt.plot([1.], [0.], marker='o', markersize=1, color="black")

plt.xlabel(r'$R/R_0$', fontsize="12")
plt.ylabel(r'$Z/R_0$',  fontsize="12")

plt.tight_layout()

plt.show()    
