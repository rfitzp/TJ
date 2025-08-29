# py.py

# Plots y component of perfect-wall ideal eigenfunction in R, Z plane.
# User prompted for solution number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('seismic')

ncont = 160

fn   = 'Equilibrium.nc'
ds   = nc.Dataset(fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']

Rw   = ds['Rwall']
Zw   = ds['Zwall']

epsa = para[0]
scale = 1.5*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'Vertical.nc'
ds1   = nc.Dataset(fn1)
psi_r = ds1['py_ideal_eig_r']

fna   = 'Equilibrium1.nc'
dsa   = nc.Dataset(fna)
Ra    = dsa['R']
Za    = dsa['Z']

Rwa   = dsa['Rwall']
Zwa   = dsa['Zwall']

RRa = np.asarray(Ra);
ZZa = np.asarray(Za);

fn1a   = 'Vertical1.nc'
ds1a   = nc.Dataset(fn1a)
psi_ra = ds1a['py_ideal_eig_r']

fig = plt.figure (figsize = (12.5, 6.0))
plt.rc ('xtick', labelsize=15) 
plt.rc ('ytick', labelsize=15) 

k = 0

Pr  = psi_r[k,:,:]
PPr = np.asarray (Pr)

plt.subplot (1, 2, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'black', linewidth = 0.5, linestyle = 'solid')

plt.plot (Rw, Zw, color = 'black', linewidth = 4, linestyle = 'solid')  

cp = plt.contourf (RR, ZZ, PPr, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "15")
plt.ylabel (r'$Z/R_0$', fontsize = "15")

plt.subplot(1, 2, 2)
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

Pra  = psi_ra[k,:,:]
PPra = np.asarray (Pra)

plt.plot (Ra[nf-1], Za[nf-1], color = 'black', linewidth = 0.5, linestyle = 'solid')    

plt.plot (Rwa, Zwa, color = 'black', linewidth = 4, linestyle = 'solid')  

ci = plt.contourf (RRa, ZZa, PPra, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color="black")

plt.xlabel (r'$R/R_0$', fontsize = "15")
plt.ylabel (r'$Z/R_0$', fontsize = "15")

plt.tight_layout()

#plt.show()    
plt.savefig ("Fig5.png")
