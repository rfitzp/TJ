# IOX.py

# Plots IO and IX versus k

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Island/Island.nc'
ds = nc.Dataset(fn)
k  = np.asarray(ds['kk'])
jo = np.asarray(ds['IO'])
jx = np.asarray(ds['IX'])
np = jo.shape[1]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Island Code: IO and IX')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (k[0], k[-1])

plt.plot (k, jo[:,0],  color = 'red',     linewidth = 2, linestyle = 'solid', label = r"$I_O(00)$")
plt.plot (k, jo[:,10], color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$I_O(10)$")
plt.plot (k, jo[:,20], color = 'blue',    linewidth = 2, linestyle = 'solid', label = r"$I_O(20)$")
plt.plot (k, jo[:,30], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = r"$I_O(30)$")
plt.plot (k, jo[:,40], color = 'magenta', linewidth = 2, linestyle = 'solid', label = r"$I_O(40)$")
plt.plot (k, jo[:,50], color = 'cyan',    linewidth = 2, linestyle = 'solid', label = r"$I_O(50)$")
plt.plot (k, jo[:,60], color = 'orange',  linewidth = 2, linestyle = 'solid', label = r"$I_O(60)$")
plt.plot (k, jo[:,70], color = 'brown',   linewidth = 2, linestyle = 'solid', label = r"$I_O(70)$")

plt.axhline ( 0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

plt.xlim (k[0], k[-1])

plt.plot (k, jx[:,0],  color = 'red',     linewidth = 2, linestyle = 'solid', label = r"$I_X(00)$")
plt.plot (k, jx[:,10], color = 'green',   linewidth = 2, linestyle = 'solid', label = r"$I_X(10)$")
plt.plot (k, jx[:,20], color = 'blue',    linewidth = 2, linestyle = 'solid', label = r"$I_X(20)$")
plt.plot (k, jx[:,30], color = 'yellow',  linewidth = 2, linestyle = 'solid', label = r"$I_X(30)$")
plt.plot (k, jx[:,40], color = 'magenta', linewidth = 2, linestyle = 'solid', label = r"$I_X(40)$")
plt.plot (k, jx[:,50], color = 'cyan',    linewidth = 2, linestyle = 'solid', label = r"$I_X(50)$")
plt.plot (k, jx[:,60], color = 'orange',  linewidth = 2, linestyle = 'solid', label = r"$I_X(60)$")
plt.plot (k, jx[:,70], color = 'brown',   linewidth = 2, linestyle = 'solid', label = r"$I_X(70)$")

plt.axhline ( 0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ( 1., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$k$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
