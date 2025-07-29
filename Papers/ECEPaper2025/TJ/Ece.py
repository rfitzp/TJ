# Ece.py

# Plots ECE convolution function parameters on tilted central chord
import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ1.nc'
ds   = nc.Dataset(fn)
R    = ds['R_eq']
Rres = ds['R_res']
sigO = ds['sigma_1^O']
sigX = ds['sigma_2^X']
delO = ds['Delta_1^O']
delX = ds['Delta_2^X']
tauO = ds['tau_1^O']
tauX = ds['tau_2^X']

fn1   = 'TJ2.nc'
ds1   = nc.Dataset(fn1)
Rres1 = ds1['R_res']

fig = plt.figure (figsize = (6.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: ECE Convolution Properties along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (3, 1, 1)

plt.xlim (1., R[-1])

plt.axvline (Rres[0],  color = 'black',  linewidth = 1.5, linestyle = 'dashed')
plt.axvline (Rres1[0], color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, sigO, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = '$\sigma_1^{\,(O)}/R_0$')
plt.plot    (R, sigX, color = 'red',   linewidth = 2,   linestyle = 'solid', label = '$\sigma_2^{\,(X)}/R_0$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R/R_0$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim (1., R[-1])

plt.axvline (Rres[0],  color = 'black',  linewidth = 1.5, linestyle = 'dashed')
plt.axvline (Rres1[0], color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, delO, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = '$\Delta_1^{\,(O)}/R_0$')
plt.plot    (R, delX, color = 'red',   linewidth = 2,   linestyle = 'solid', label = '$\Delta_2^{\,(X)}/R_0$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R/R_0$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (1., R[-1])

plt.axvline (Rres[0],  color = 'black',  linewidth = 1.5, linestyle = 'dashed')
plt.axvline (Rres1[0], color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, tauO, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r'$\tau_{\infty\,1}^{\,(O)}$')
plt.plot    (R, tauX, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r'$\tau_{\infty\,2}^{\,(X)}$')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R/R_0$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Fig15.pdf")
