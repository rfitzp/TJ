# Psi.py

# Plots Psi component of unreconnected eigenfunction associated with given rational surface in R, Z plane.
# User prompted for rational surface number and scaling exponent

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
scale = 1.2*epsa

RR = np.asarray(R);
ZZ = np.asarray(Z);
rr = np.asarray(r);
nf = RR.shape[0]
nt = RR.shape[1]                

fn1   = 'TJ1.nc'
ds1   = nc.Dataset(fn1)
rres  = ds1['r_res']
psi_r = ds1['Psi_unrc_eig_r']
psi_i = ds1['Psi_unrc_eig_i']

fig = plt.figure (figsize = (8, 7.5))
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

"""
m = input ("rational surface number ? ")
k = int(m) - 1
"""

k = 0

"""
scl = input ("scaling exponent ? ")
"""

scl = 1

Pr  = -psi_r[k,:,:]
Pi  = -psi_i[k,:,:]

PPr = np.asarray (Pr)
PPi = np.asarray (Pi)

prmin = np.amin (PPr)
prmax = np.amax (PPr)
pimin = np.amin (PPi)
pimax = np.amax (PPi)

for i in range (np.size(PPr,0)):
    for j in range (np.size(PPr,1)):
        val = PPr[i][j]
        if (val > 0.):
            val1 = float(val /prmax)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /prmin))**float(scl)
        else:
            val1 = 0.
        PPr[i][j] = val1

for i in range (np.size(PPi,0)):
    for j in range (np.size(PPi,1)):
        val = PPi[i][j]
        if (val > 0.):
            val1 = float(val /pimax)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /pimin))**float(scl)
        else:
            val1 = 0.
        PPi[i][j] = val1

Pr1  = psi_r[1,:,:]
Pi1  = psi_i[1,:,:]

PPr1 = np.asarray (Pr1)
PPi1 = np.asarray (Pi1)

prmin1 = np.amin (PPr1)
prmax1 = np.amax (PPr1)
pimin1 = np.amin (PPi1)
pimax1 = np.amax (PPi1)

for i in range (np.size(PPr1,0)):
    for j in range (np.size(PPr1,1)):
        val = PPr1[i][j]
        if (val > 0.):
            val1 = float(val /prmax1)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /prmin1))**float(scl)
        else:
            val1 = 0.
        PPr1[i][j] = val1

for i in range (np.size(PPi1,0)):
    for j in range (np.size(PPi1,1)):
        val = PPi1[i][j]
        if (val > 0.):
            val1 = float(val /pimax1)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /pimin1))**float(scl)
        else:
            val1 = 0.
        PPi1[i][j] = val1
        
plt.subplot (1, 1, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-0.95*scale, 0.95*scale)
plt.axis('off')

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

cp = plt.contourf (RR, ZZ, PPr, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout()

#plt.show()    
plt.savefig("Psi2.jpeg")
#plt.savefig("Psi1a.png")
