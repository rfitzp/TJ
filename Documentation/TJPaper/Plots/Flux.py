# Flux.py

# Plots r, theta coordinate system in R, Z plane.
# Also shows rational surfaces and RMP coils.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'Equilibrium.nc'
ds   = nc.Dataset (fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']

epsa = para[0]

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'TJ.nc'
ds1   = nc.Dataset (fn1)
rres  = ds1['rres']
Rcoil = np.asarray(ds1['Rcoil'])
Zcoil = np.asarray(ds1['Zcoil'])
Icoil = np.asarray(ds1['Icoil'])

fig = plt.figure (figsize = (8.0, 7.5))
fig.canvas.manager.set_window_title (r'TJ Code: r, theta Coordinate System')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17) 

plt.subplot(1, 1, 1)

scale = 1.5*epsa
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

#for n in range (0, nf):
#    plt.plot (R[n], Z[n],   color = 'blue', linewidth = 0.1, linestyle = 'solid')
#plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

#for n in range (0, nt-1):
#    plt.plot (R[:,n], Z[:,n], color = 'green', linewidth = 0.1, linestyle = 'solid')

for n in range (0, nf, 10):
    plt.plot (R[n], Z[n], color = 'blue', linewidth = 0.5, linestyle = 'solid')

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

for n in range (0, nt-1, 10):
    plt.plot (R[:,n], Z[:,n], color = 'green', linewidth = 0.5, linestyle = 'solid')

plt.contour (RR, ZZ, rr, rres, colors = 'red', linewidths = 1.)    

plt.plot ([1.], [0.], marker = 'o', markersize = 2, color = "black")

for R, Z, I in zip (Rcoil, Zcoil, Icoil):
    if I > 0.:
        plt.plot ([R], [Z], marker = 'o', markersize = 7, color = "blue")
    if I < 0.:
        plt.plot ([R], [Z], marker = 'o', markersize = 7, color = "red")

plt.xlabel (r'$R/R_0$', fontsize="20")
plt.ylabel (r'$Z/R_0$', fontsize="20")

plt.tight_layout ()

plt.savefig ("Figure2.pdf")
#plt.show ()
