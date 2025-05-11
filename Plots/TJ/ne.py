# ne.py

# Plots n_e component of unreconnected eigenfunction associated with given rational surface in R, Z plane.
# User prompted for rational surface number and toroidal gridpoint.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
ReBu = plt.get_cmap ('rainbow') #plt.get_cmap ('seismic')

ncont = 360

fn   = '../../Outputs/Equilibrium/Equilibrium.nc'
ds   = nc.Dataset (fn)
para = ds['para']
R    = ds['R']
Z    = ds['Z']
r    = ds['rr']
t    = ds['theta']

epsa = para[0]
scale = 1.5*epsa

RR = np.asarray (R);
ZZ = np.asarray (Z);
rr = np.asarray (r);
nf = RR.shape[0]
nt = RR.shape[1]

fn1   = '../../Outputs/TJ/TJ.nc'
ds1   = nc.Dataset (fn1)
rres  = ds1['r_res']
psi_r = ds1['ne']
pn    = psi_r.shape[3]  

fig = plt.figure (figsize = (8.5, 8.))
fig.canvas.manager.set_window_title (r'TJ Code: n_e(R, Z)')
plt.rc ('xtick', labelsize=12) 
plt.rc ('ytick', labelsize=12) 

nres = len (rres)
m = input ("rational surface number (%d .. %d) ? " % (1, nres))
k = int (m) - 1
nn = input ("toroidal gridpoint (%d .. %d) ? " % (1, pn))
n1 = int (nn) - 1

Pr0 = np.asarray (psi_r[k,:,:,n1])

plt.subplot (1, 1, 1)
plt.xlim (1.-scale, 1.+scale)
plt.ylim (-scale, scale)

plt.plot (R[nf-1], Z[nf-1], color = 'blue', linewidth = 0.5, linestyle = 'solid')    

plt.contour (RR, ZZ, Pr0, levels = 360, colors = 'white', linewidths = 0.2)

cp = plt.contourf (RR, ZZ, Pr0, ncont, cmap = ReBu)

plt.contour (RR, ZZ, rr, rres, colors = 'black', linewidths = 0.5)    

plt.plot ([1.], [0.], marker = 'o', markersize = 1, color = "black")

plt.xlabel (r'$R/R_0$', fontsize = "12")
plt.ylabel (r'$Z/R_0$', fontsize = "12")

plt.tight_layout ()

plt.show ()    
#plt.savefig("dTe.png")
