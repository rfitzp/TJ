# deltaW.py

# Plots delta_W versus solution number

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator

fn = '../../Outputs/TJ/TJ.nc'
ds = nc.Dataset(fn)
w  = ds['delta_W']

J = len(w)

jj = np.linspace (0, J-1, J-1, endpoint = False)

ww = np.asarray(w)

www = []
for j in range (J-1):
    www.append(ww[j])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: delta W Values')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (-1, J-1)

plt.plot    (jj, www, color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10)
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

ax = fig.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel (r'$j$',                fontsize = "15")
plt.ylabel (r'$\delta \hat{W}_j$', fontsize = "15")


plt.tight_layout ()

plt.show ()    
