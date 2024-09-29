# PsiZr.py

# Plots Psi, Z components of ideal RMP response in R, Z plane.

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
psi_r = ds1['Psi_rmp_eig_r']
psi_i = ds1['Psi_rmp_eig_i']
z_r   = ds1['Z_rmp_eig_r']
z_i   = ds1['Z_rmp_eig_i']

fig = plt.figure (figsize = (8.5, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Psi_rmp(R, Z), Z_rmp(R, Z)')
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 


scl = input ("scaling exponent ? ")

Pr  = psi_r[:,:]
Pi  = psi_i[:,:]
Zr  = z_r[:,:]
Zi  = z_i[:,:]

PPr = np.asarray (Pr)
PPi = np.asarray (Pi)
ZZr = np.asarray (Zr)
ZZi = np.asarray (Zi)

prmin = np.amin (PPr)
prmax = np.amax (PPr)
pimin = np.amin (PPi)
pimax = np.amax (PPi)
zrmin = np.amin (ZZr)
zrmax = np.amax (ZZr)
zimin = np.amin (ZZi)
zimax = np.amax (ZZi)

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

for i in range (np.size(ZZr,0)):
    for j in range (np.size(ZZr,1)):
        val = ZZr[i][j]
        if (val > 0.):
            val1 = float(val /zrmax)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /zimin))**float(scl)
        else:
            val1 = 0.
        ZZr[i][j] = val1

for i in range (np.size(ZZi,0)):
    for j in range (np.size(ZZi,1)):
        val = ZZi[i][j]
        if (val > 0.):
            val1 = float(val /zimax)**float(scl)
        elif (val < 0.):
            val1 = - (float(val /zimin))**float(scl)
        else:
            val1 = 0.
        ZZi[i][j] = val1   
        
plt.subplot (2, 2, 1)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors='black', linewidths = 0.5)    

plt.contourf (RR, ZZ, Pr, ncont, cmap = ReBu)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot (2, 2, 2)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

plt.contourf (RR, ZZ, PPi, ncont, cmap=ReBu)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot (2, 2, 3)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

plt.contourf (RR, ZZ, ZZr, ncont, cmap = ReBu)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot (2, 2, 4)

plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

plt.contourf (RR, ZZ, ZZi, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout()

plt.show()    
