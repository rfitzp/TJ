import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn     = 'TearX.nc'
ds     = nc.Dataset(fn)
r      = np.asarray(ds['r'])
omegae = ds['omegae']
omegai = ds['omegai']
omegaE = ds['omegaE']

oe = np.asarray(omegae) + np.asarray(omegaE)

for i in range(len(oe)-1):
    if oe[i]*oe[i+1] < 0.:
        x = (r[i]*oe[i+1]-r[i+1]*oe[i])/(oe[i+1]-oe[i])

print ("r_0 = %11.4e" % x)        

fontsize = 17

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0.8, 1.)
 
plt.plot (r, omegae, color = 'red',   linewidth = 2,    linestyle = 'solid', label = r"$\omega_{\ast\,e}$")
plt.plot (r, omegai, color = 'blue',  linewidth = 2,   linestyle = 'solid',  label = r"$\omega_{\ast\,i}$")
plt.plot (r, omegaE, color = 'green', linewidth = 2,   linestyle = 'solid',  label = r"$\omega_E$")
plt.plot (r, oe,     color = 'black', linewidth = 2,   linestyle = 'solid',  label = r"$\omega_{\perp\,e}$")
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (x,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))    
plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure8.pdf")
