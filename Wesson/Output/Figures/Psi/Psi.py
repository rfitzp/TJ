# Psi.py

# Plots tearing eigenfunction versus radius
# User promted for rational surface number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = 'Tear.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
psi = ds['psi']
rs  = ds['rres']

fn1  = 'Tear1.nc'
ds1  = nc.Dataset(fn1)
r1   = ds1['r']
psi1 = ds1['psi']
rs1  = ds1['rres']

p2 = np.asarray(psi[:,0])/2.
p3 = np.asarray(psi[:,1])/3.

p2a = np.asarray(psi1[:,0])/2.
p3a = np.asarray(psi1[:,1])/3.

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, p2, color = 'black', linewidth = 3,   linestyle = 'solid',  label = "$m=2/n=1$")
plt.plot    (r, p3, color = 'black', linewidth = 3,   linestyle = 'dashed', label = "$m=3/n=1$")
plt.axhline (0.,          color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs[0],       color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axvline (rs[1],       color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axhline (1.,          color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$',    fontsize = "20")
plt.ylabel (r'$\hat{\psi}$', fontsize = "20")
plt.title (r"$\bar{b} = 10$", fontsize = 15)
plt.legend (fontsize = 15)

plt.subplot (1, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r, p2a, color = 'black', linewidth = 3,   linestyle = 'solid',  label = "$m=2/n=1$")
plt.plot    (r, p3a, color = 'black', linewidth = 3,   linestyle = 'dashed', label = "$m=3/n=1$")
plt.axhline (0.,          color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (rs[0],       color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axvline (rs[1],       color = 'black', linewidth = 2.0, linestyle = 'dotted')
plt.axhline (1.,          color = 'black', linewidth = 2.0, linestyle = 'dotted')

plt.xlabel (r'$\bar{r}$',    fontsize = "20")
plt.ylabel (r'$\hat{\psi}$', fontsize = "20")
plt.title (r"$\bar{b} = 1.05$", fontsize = 15)
plt.legend (fontsize = 15)

plt.tight_layout ()

plt.show ()    

#plt.savefig ("Figure11_1.pdf")
