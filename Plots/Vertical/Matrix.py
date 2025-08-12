# Matrix.py

# Plots A_m^m', B_m^m', C_m^m', and D_m^m' coupling matrices versus radius.
# User is prompted for values of m and m'.
# Also plots symmetry tests for matrices versus r.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = '../../Outputs/Vertical/Vertical.nc'
ds     = nc.Dataset(fn)
r      = ds['r']
mpol   = ds['mpol']
lmmpr  = ds['Ammp_r']
lmmpi  = ds['Ammp_i']
mmmpr  = ds['Bmmp_r']
mmmpi  = ds['Bmmp_i']
nmmpr  = ds['Cmmp_r']
nmmpi  = ds['Cmmp_i']
pmmpr  = ds['Dmmp_r']
pmmpi  = ds['Dmmp_i']
ltest  = ds['Atest']
mntest = ds['BCtest']
ptest  = ds['Dtest']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Coupling Matrices')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

m   = input ("m  (%d .. %d) ? " % (mpol[0], mpol[-1]))
mp  = input ("m' (%d .. %d) ? " % (mpol[0], mpol[-1]))
j   = int(m)  - mpol[0]
jp  = int(mp) - mpol[0]
if j < 0 or jp < 0:
    print ("m m' out of range")
    quit()
if j > lmmpr.shape[1] or jp > lmmpr.shape[1]:
    print ("m m' out of range")
    quit()    

plt.subplot (3, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, lmmpr[:,j,jp], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r, lmmpi[:,j,jp], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,               color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r"$A_m^{m'}$", fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (0., 1.)

plt.plot    (r,  mmmpr[:,j,jp], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r,  mmmpi[:,j,jp], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r"$B_m^{m'}$", fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 1.)

plt.plot    (r,  nmmpr[:,j,jp], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r,  nmmpi[:,j,jp], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r"$C_m^{m'}$", fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 1.)

plt.plot    (r,  pmmpr[:,j,jp], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r,  pmmpi[:,j,jp], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',  fontsize = "15")
plt.ylabel (r"$D_m^{m'}$", fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim(0., 1.)
plt.ylim(-1.e-14, 1.e-14)

plt.plot    (r,  ltest[:,j,jp],  color = 'blue',   linewidth = 2,   linestyle = 'solid')
plt.plot    (r,  mntest[:,j,jp], color = 'green',  linewidth = 2,   linestyle = 'solid')
plt.plot    (r,  ptest[:,j,jp],  color = 'yellow', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                 color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',      fontsize = "15")
plt.ylabel (r"Symmetry tests", fontsize = "15")

plt.tight_layout ()

plt.show ()    
