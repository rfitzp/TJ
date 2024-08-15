# h.py

# Plots adaptive ode integration step-length, h, and truncation error, err, versus radius, r.
# Locations of rational surfaces are shown.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ.nc'
ds   = nc.Dataset(fn)
r    = ds['r_grid']
h    = ds['h_ode']
err  = ds['err_ode']
rres = ds['r_res']

rr   = np.asarray(r)
hh   = np.asarray(h)
eerr = np.asarray(err)

hhh   = np.log10(hh);
eeerr = np.log10(eerr);

fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Cash-Karp RK4/RK5 Adaptive Integration')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.xlim (0., 1.)

for rx in rres:
    plt.axvline (rx, color = 'red', linewidth = 1.5, linestyle = 'dashed')

plt.plot (rr, hhh, color = 'blue', linewidth = 1, linestyle = 'solid')

plt.xlabel (r'$\hat{r}$',      fontsize = "15")
plt.ylabel (r'$\log_{10}(h)$', fontsize = "15")

plt.subplot(2, 1, 2)

plt.xlim(0., 1.)

for rx in rres:
    plt.axvline (rx, color='red', linewidth = 1.5, linestyle = 'dashed')

plt.plot    (rr,  eeerr, color = 'blue', linewidth = 1,    linestyle = 'solid')
plt.axhline (-12,        color = 'red',  linewidth = 0.5,  linestyle = 'dotted')
plt.axhline (-13,        color = 'red',  linewidth = 0.5,  linestyle = 'dotted')
plt.axhline (-14,        color = 'red',  linewidth = 0.5,  linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',             fontsize = "15")
plt.ylabel (r'$\log_{10}(\epsilon)$', fontsize = "15")

plt.tight_layout ()

plt.show ()    
