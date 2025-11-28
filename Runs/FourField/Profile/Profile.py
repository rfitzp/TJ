import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TearX.nc'
ds   = nc.Dataset(fn)
r    = ds['r']
q    = ds['q']
s    = ds['s']
ne   = ds['ne']
Te   = ds['Te']
Ti   = ds['Ti']

ss = np.asarray(s)
s1 = 1. + ss
ls = np.log(s1)

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 2, 1)

plt.xlim (0.8, 1.)
plt.ylim (2.,  9.)
 
plt.plot (r, q, color = 'blue',  linewidth = 2,   linestyle = 'solid')

plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.ylabel (r'$q$',       fontsize = fontsize)

plt.subplot (2, 2, 2)

plt.xlim (0.8, 1.)
 
plt.plot    (r, ls, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',  fontsize = fontsize)
plt.ylabel (r'$\ln(1+s)$', fontsize = fontsize)

plt.subplot (2, 2, 3)

plt.xlim (0.8, 1.)
plt.ylim (0., 1.5e20)

plt.plot (r, ne, color = 'blue',  linewidth = 2,   linestyle = 'solid')

plt.xlabel (r'$\hat{r}$',              fontsize = fontsize)
plt.ylabel (r'$n_e(10^{19}\,m^{-3})$', fontsize = fontsize)

plt.subplot (2, 2, 4)

plt.xlim (0.8, 1.)
plt.ylim (0., 7.e3)
 
plt.plot (r, Te, color = 'red',  linewidth = 2,   linestyle = 'solid', label = "$T_e\,\,(eV)$")
plt.plot (r, Ti, color = 'blue', linewidth = 2,   linestyle = 'solid', label = "$T_i\,\,(eV)$")

plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))    
plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure7.pdf")
