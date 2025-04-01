# dTe.py

# Plots delta T_e component of unreconnected eigenfunction associated with given rational surface in R, Z plane.
# User prompted for rational surface number and scaling exponent

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
psi_r = ds1['dTe_cos']
psi_i = ds1['dTe_sin']

fig = plt.figure (figsize = (12.5, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: delta T_e(R, Z)')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m = input ("rational surface number ? ")
k = int(m) - 1

scl = input ("scaling exponent ? ")

Pr  = psi_r[k,:,:]
Pi  = psi_i[k,:,:]

PPr = np.asarray (Pr)
PPi = np.asarray (Pi)

prmin = np.amin (PPr)
prmax = np.amax (PPr)
pimin = np.amin (PPi)
pimax = np.amax (PPi)

print ("delta T_e_cos = (%10.3e, %10.3e)" % (prmin, prmax))
print ("delta T_e_sin = (%10.3e, %10.3e)" % (pimin, pimax))

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

plt.subplot (1, 2, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

cp = plt.contourf (RR, ZZ, PPr, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.subplot(1, 2, 2)
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

ci = plt.contourf (RR, ZZ, PPi, ncont, cmap = ReBu)

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color="black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout()

plt.show()    
