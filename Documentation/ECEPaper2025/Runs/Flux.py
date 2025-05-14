# Flux.py

# Plots r, theta coordinate system in R, Z plane.
# Also shows rational surfaces, wall, and RMP coils.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'Equilibrium1.nc'
ds   = nc.Dataset (fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']
Req  = ds['R_eq']
Zeq  = ds['Z_eq']

epsa = para[0]

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'TJ1.nc'
ds1   = nc.Dataset (fn1)
rres  = np.asarray(ds1['r_res'])

rrex = rres[0]

fn2   = 'TJ2.nc'
ds2   = nc.Dataset (fn2)
rre2  = np.asarray(ds2['r_res'])

fig = plt.figure (figsize = (8.0, 7.5))
fig.canvas.manager.set_window_title (r'TJ Code: r, theta Coordinate System')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17) 

plt.subplot(1, 1, 1)

scale = 1.3*epsa
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

for n in range (0, nf, 10):
    plt.plot (R[n], Z[n],   color = 'blue', linewidth = 0.1, linestyle = 'solid')
plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

for n in range (0, nt-1, 10):
    plt.plot (R[:,n], Z[:,n], color = 'green', linewidth = 0.1, linestyle = 'solid')

for n in range (0, nf, 10):
    plt.plot (R[n], Z[n], color = 'blue', linewidth = 0.5, linestyle = 'solid')

plt.plot (R[-1], Z[-1], color = 'blue', linewidth = 1.0, linestyle = 'solid')    

for n in range (0, nt-1, 10):
    plt.plot (R[:,n], Z[:,n], color = 'green', linewidth = 0.5, linestyle = 'solid')

plt.contour (RR, ZZ, rr, levels=[rre2[0], rres[0]], colors = 'red', linewidths = 1.)

plt.plot (Req, Zeq, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.plot ([1.], [0.], marker = 'o', markersize = 2, color = "black")


plt.xlabel (r'$R/R_0$', fontsize = "20")
plt.ylabel (r'$Z/R_0$', fontsize = "20")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Fig1.pdf")
