# ece.py

# Plots delta T_e and delta n_e components of unreconnected eigenfunction associated
# with given rational surface versus r at given angular gridpoint.
# User prompted for rational surface number and gridpoint.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Equilibrium/Equilibrium.nc'
ds  = nc.Dataset(fn)
rr = ds['rr']

fn1 = '../../Outputs/TJ/TJ.nc'
ds1 = nc.Dataset(fn1)
dnc = ds1['dne_cos']
dns = ds1['dne_sin']
dTc = ds1['dTe_cos']
dTs = ds1['dTe_sin']
rres = ds1['r_res']

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m = input ("rational surface number ? ")
k = int(m) - 1

j = dnc.shape[2]
h = input ("grid point (0, %4d) ? " % (j-1))
g = int(h)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r"TJ Code: Components of delta T_e and delta n_e at Angular Gridpoint %d" % g)
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (rr[:,0], dTc[k,:,g], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (rr[:,0], dTs[k,:,g], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,                  color = 'black',  linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r'$\delta T_e(eV)$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 1.)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (rr[:,0], dnc[k,:,g], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (rr[:,0], dns[k,:,g], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,                  color = 'black',  linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r'$\delta n_e(m^{-3})$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
#plt.savefig("ece.png")
