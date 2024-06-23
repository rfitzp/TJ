# PsiZ.py

# Plots Psi component of unreconnected eigenfunction associated with given rational surface in R, Z plane.
# User prompted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap('seismic')

ncont = 160

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
para = ds['para']
R  = ds['R']
Z  = ds['Z']
r  = ds['rr']
t  = ds['theta']

epsa = para[0]
scale = 1.5*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'TJ.nc'
ds1   = nc.Dataset(fn1)
rres  = ds1['rres']
psi_r = ds1['Psi_unrc_eig_r']
psi_i = ds1['Psi_unrc_eig_i']
z_r   = ds1['Z_unrc_eig_r']
z_i   = ds1['Z_unrc_eig_i']

fig = plt.figure(figsize=(12.0, 6.0))
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 

m = input ("rational surface number ? ")
k = int(m) - 1

Pr  = psi_r[k,:,:]
PPr = np.asarray(Pr)
Pi  = psi_i[k,:,:]
PPi = np.asarray(Pi)
Zr  = z_r[k,:,:]
ZPr = np.asarray(Zr)
Zi  = z_i[k,:,:]
ZPi = np.asarray(Zi)

plt.subplot(1, 2, 1)
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

plt.plot(R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour(RR, ZZ, rr, rres, colors='black', linewidths = 0.5)    

plt.contourf(RR, ZZ, Pr, ncont, cmap=ReBu)    

plt.plot([1.], [0.], marker='o', markersize=1, color="black")

plt.xlabel(r'$R/R_0$', fontsize="12")
plt.ylabel(r'$Z$',  fontsize="12")

plt.subplot(1, 2, 2)
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

plt.plot(R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour(RR, ZZ, rr, rres, colors='black', linewidths = 0.5)    

plt.contourf(RR, ZZ, Pi, ncont, cmap=ReBu)    

plt.plot([1.], [0.], marker='o', markersize=1, color="black")

plt.xlabel(r'$R/R_0$', fontsize="12")
plt.ylabel(r'$Z$',  fontsize="12")

plt.tight_layout()

plt.show()    
