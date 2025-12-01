# c_beta comparison

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn  = 'FourField020.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = - np.asarray(ds['g_i'])
ds  = nc.Dataset(fn)
d4r = np.asarray(ds['Delta4_r'])
d4i = - np.asarray(ds['Delta4_i'])

fn1  = 'FourField010.nc'
ds1  = nc.Dataset(fn1)
d4r1 = np.asarray(ds1['Delta4_r'])
d4i1 = - np.asarray(ds1['Delta4_i'])

fn2  = 'FourField001.nc'
ds2  = nc.Dataset(fn2)
d4r2 = np.asarray(ds2['Delta4_r'])
d4i2 = - np.asarray(ds2['Delta4_i'])
d3r  = np.asarray(ds2['Delta3_r'])
d3i  = - np.asarray(ds2['Delta3_i'])

fn3  = 'FourField005.nc'
ds3  = nc.Dataset(fn3)
d4r3 = np.asarray(ds3['Delta4_r'])
d4i3 = - np.asarray(ds3['Delta4_i'])

dr100 = d4r  - d3r
di100 = d4i  - d3i
dr010 = d4r1 - d3r
di010 = d4i1 - d3i
dr001 = d4r2 - d3r
di001 = d4i2 - d3i
dr005 = d4r3 - d3r
di005 = d4i3 - d3i

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, dr100, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.2$")
plt.plot (gi, dr010, color = 'green', linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.1$")
plt.plot (gi, dr005, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.05$")
plt.plot (gi, dr001, color = 'black', linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.01$")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$',                                      fontsize = fontsize)
plt.ylabel (r'${\rm Re} (\hat{\Delta}_4-\hat{\Delta}_3)$', fontsize = fontsize)
plt.legend (fontsize = 15);

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, di100, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.2$")
plt.plot (gi, di010, color = 'green', linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.1$")
plt.plot (gi, di005, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.05$")
plt.plot (gi, di001, color = 'black', linewidth = 2, linestyle = 'solid', label = r"$c_\beta=0.01$")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$',                                      fontsize = fontsize)
plt.ylabel (r'${\rm Im} (\hat{\Delta}_4-\hat{\Delta}_3)$', fontsize = fontsize)
plt.legend (fontsize = 15);

plt.tight_layout ()

#plt.show () 
plt.savefig ("Figure5.pdf")

