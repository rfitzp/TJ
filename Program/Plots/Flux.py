# Flux.py

# Plots equilibrium magnetic flux-surfaces in R, Z plane.
# User prompted for scale (x-limits: 1.-scale, 1.+scale; z-limits: -scale, scale)

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds = nc.Dataset(fn)
R  = ds['R']
Z  = ds['Z']

RR = np.asarray(R)
nf = RR.shape[0]
nt = RR.shape[1]                

fig = plt.figure(figsize=(8.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)

scale = float (input ("scale ? "))
plt.xlim(1.-scale, 1.+scale)
plt.ylim(-scale, scale)

for n in range (nf-1):
    plt.plot(R[n], Z[n], color = 'blue', linewidth = 0.5, linestyle = 'solid')
plt.plot(R[nf-1],  Z[nf-1],  color = 'red', linewidth = 1.0, linestyle = 'solid')

for n in range (0,nt-1,5):
    plt.plot(R[:,n], Z[:,n], color = 'green', linewidth = 0.5, linestyle = 'solid')

plt.plot([1.], [0.], marker='o', markersize=2, color="red")

plt.xlabel(r'$R/R_0$', fontsize="20")
plt.ylabel(r'$Z/Z_0$',  fontsize="20")

plt.tight_layout()

plt.show()    
