# Flux.py

# Plots equilibrium magnetic flux-surfaces in R, Z plane.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
para = ds['para']
R  = ds['R']
Z  = ds['Z']
r  = ds['rr']
t  = ds['theta']

epsa = para[0]

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1 = 'TJ.nc'
ds1 = nc.Dataset(fn1)
rres = ds1['rres']

fig = plt.figure(figsize=(8.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)

scale = 1.5*epsa
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

for n in range (0, nf, 10):
    plt.plot(R[n], Z[n], color = 'blue', linewidth = 0.5, linestyle = 'solid')
plt.plot(R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

for n in range (0,nt-1,5):
    plt.plot(R[:,n], Z[:,n], color = 'green', linewidth = 0.5, linestyle = 'solid')

plt.contour(RR, ZZ, rr, rres, colors='red', linewidths = 1.)    

plt.plot([1.], [0.], marker='o', markersize=2, color="red")

plt.xlabel(r'$R/R_0$', fontsize="20")
plt.ylabel(r'$Z$',  fontsize="20")

plt.tight_layout()

plt.show()    
