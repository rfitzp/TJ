# km.py

# Plots k_m versus radius.
# User is prompted for values of m.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = '../../Outputs/TJ/TJ.nc'
ds     = nc.Dataset(fn)
r      = ds['r']
mpol   = ds['mpol']
rres   = ds['r_res']
km     = ds['km']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: k_m')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

m   = input ("m  (%d .. %d)? " % (mpol[0], mpol[-1]))
j   = int(m)  - mpol[0]
if j < 0:
    print ("m out of range")
    quit()
if j > km.shape[1]:
    print ("m out of range")
    quit()    

plt.subplot (1, 1, 1)

plt.xlim (0., 1.)

plt.plot    (r, km[:,j], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,         color = 'black', linewidth = 1.5, linestyle = 'dotted')

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r"$k_m$",     fontsize = "15")

plt.tight_layout ()

plt.show ()    
