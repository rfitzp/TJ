# Fscan.py

# Plots real and imaginary parts of Delta3 and Delta4 versus imagninary part of g

# Qe = 1, Qi = -1, D = 0.5, Pphi = 1, Pperp = 1, JK = 0

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Run1.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = -np.asarray(ds['g_i'])
d3r = ds['Delta3_r']
d3i = ds['Delta3_i']
d4r = ds['Delta4_r']
d4i = ds['Delta4_i']

fn1 = 'Run2.nc'
ds1  = nc.Dataset(fn1)
In1  = ds1['InputParameters']
gi1  = -np.asarray(ds1['g_i'])
d3r1 = ds1['Delta3_r']
d3i1 = ds1['Delta3_i']
d4r1 = ds1['Delta4_r']
d4i1 = ds1['Delta4_i']

fn2 = 'Run3.nc'
ds2  = nc.Dataset(fn2)
In2  = ds2['InputParameters']
gi2  = -np.asarray(ds2['g_i'])
d3r2 = ds2['Delta3_r']
d3i2 = ds2['Delta3_i']
d4r2 = ds2['Delta4_r']
d4i2 = ds2['Delta4_i']


fig = plt.figure (figsize = (10.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.plot (gi, d3r,  color = 'black',  linewidth = 2, linestyle = 'dotted',  label = r"$c_{\beta\,k} = 0.0$")
plt.plot (gi, d4r,  color = 'black',  linewidth = 2, linestyle = 'dashed',  label = r"$c_{\beta\,k} = 0.1$")
plt.plot (gi, d4r2, color = 'black',  linewidth = 2, linestyle = 'dashdot', label = r"$c_{\beta\,k} = 0.5$")
plt.plot (gi, d4r1, color = 'black',  linewidth = 2, linestyle = 'solid',   label = r"$c_{\beta\,k} = 1.0$")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_k-Q_{E\,k}$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.plot (gi, d3i,  color = 'black',  linewidth = 2, linestyle = 'dotted', label = r"$c_{\beta\,k} = 0.0$")
plt.plot (gi, d4i,  color = 'black',  linewidth = 2, linestyle = 'dashed', label = r"$c_{\beta\,k} = 0.1$")
plt.plot (gi, d4i2, color = 'black',  linewidth = 2, linestyle = 'dashdot', label = r"$c_{\beta\,k} = 0.5$")
plt.plot (gi, d4i1, color = 'black',  linewidth = 2, linestyle = 'solid',  label = r"$c_{\beta\,k} = 1.0$")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_k-Q_{E\,k}$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

#plt.show () 
plt.savefig ("Fig8_3.pdf")
