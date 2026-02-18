# Scan.py

# Plots results of qa scan

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/Kink/Kink.nc'
ds   = nc.Dataset(fn)
p    = ds['para']
qa   = ds['q_a']
ll   = ds['lambda']
gg   = ds['gamma']
bb   = ds['b_crit']

rs = p[0];

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'KINK Code: q_a Scan: m = %2d n = %2d nu = %10.3e' % (p[1], p[2], p[3]))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 1, 1)

plt.plot (qa, ll, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,  color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$q_a$',     fontsize = "15")
plt.ylabel (r'$\lambda$', fontsize = "15")

plt.subplot (3, 1, 2)

plt.plot (qa, gg, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,  color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$q_a$',            fontsize = "15")
plt.ylabel (r'$\gamma\,\tau_w$', fontsize = "15")

plt.subplot (3, 1, 3)

plt.plot (qa, bb, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (1.,  color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$q_a$',            fontsize = "15")
plt.ylabel (r'$\bar{b}_{crit}$', fontsize = "15")


plt.tight_layout ()

plt.show ()    
