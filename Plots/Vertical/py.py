# py.py

# Plots y component of perfect-wall ideal eigenfunction in R, Z plane.
# User prompted for solution number

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

Rw   = ds['Rwall']
Zw   = ds['Zwall']

epsa = para[0]
scale = 1.5*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = '../../Outputs/Vertical/Vertical.nc'
ds1   = nc.Dataset(fn1)
psi_r = ds1['py_ideal_eig_r']
psi_i = ds1['py_ideal_eig_i']
w     = ds1['pdelta_W']

J = len(w)

fig = plt.figure (figsize = (12.5, 6.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Perfect-Wall y(R, Z)')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

print ("solution number = (%d .. %d)" % (0, J-1))
m = input ("solution number ? ")
k = int(m)

Pr  = psi_r[k,:,:]
Pi  = psi_i[k,:,:]

PPr = np.asarray (Pr)
PPi = np.asarray (Pi)

plt.subplot (1, 2, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')

plt.plot (Rw, Zw, color = 'black', linewidth = 4, linestyle = 'solid')  

cp = plt.contourf (RR, ZZ, PPr, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot(1, 2, 2)
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.plot (Rw, Zw, color = 'black', linewidth = 4, linestyle = 'solid')  

ci = plt.contourf (RR, ZZ, PPi, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color="black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout()

plt.show()    
