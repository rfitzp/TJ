# Qe = 1, Qi = -1, D = 2.0, Pphi = 1, Pperp = 1, JK = 0

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

delta = 0.1

fn = 'Run1a.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = -np.asarray(ds['g_i'])
d3r = np.asarray(ds['Delta3_r'])
d3i = np.asarray(ds['Delta3_i'])
d4r = np.asarray(ds['Delta4_r'])
d4i = np.asarray(ds['Delta4_i'])

T3 = delta * d3i /((d3r + delta)*(d3r + delta) + d3i*d3i)
T4 = delta * d4i /((d4r + delta)*(d4r + delta) + d4i*d4i)

fn1 = 'Run2a.nc'
ds1  = nc.Dataset(fn1)
d4r1 = np.asarray(ds1['Delta4_r'])
d4i1 = np.asarray(ds1['Delta4_i'])

T1 = delta * d4i1 /((d4r1 + delta)*(d4r1 + delta) + d4i1*d4i1)

fn2 = 'Run3a.nc'
ds2  = nc.Dataset(fn2)
d4r2 = np.asarray(ds2['Delta4_r'])
d4i2 = np.asarray(ds2['Delta4_i'])

T2 = delta * d4i2 /((d4r2 + delta)*(d4r2 + delta) + d4i2*d4i2)

fig = plt.figure (figsize = (8.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (3, 1, 1)

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
plt.xlabel (r'$Q_k-Q_{E\,k}$',        fontsize = "15")
plt.ylabel (r'$Re(\hat{\Delta}_k)}$', fontsize = "15")
plt.legend (fontsize = "12");

plt.subplot (3, 1, 2)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.plot (gi, d3i,  color = 'black',  linewidth = 2, linestyle = 'dotted',  label = r"$c_{\beta\,k} = 0.0$")
plt.plot (gi, d4i,  color = 'black',  linewidth = 2, linestyle = 'dashed',  label = r"$c_{\beta\,k} = 0.1$")
plt.plot (gi, d4i2, color = 'black',  linewidth = 2, linestyle = 'dashdot', label = r"$c_{\beta\,k} = 0.5$")
plt.plot (gi, d4i1, color = 'black',  linewidth = 2, linestyle = 'solid',   label = r"$c_{\beta\,k} = 1.0$")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_k-Q_{E\,k}$',        fontsize = "15")
plt.ylabel (r'$Im(\hat{\Delta}_k)}$', fontsize = "15")
plt.legend (fontsize = "12");

plt.subplot (3, 1, 3)

plt.xlim (gi[0], 0)
plt.gca().invert_xaxis()

plt.plot (gi, T3, color = 'black',  linewidth = 2, linestyle = 'dotted',  label = r"$c_{\beta\,k} = 0.0$")
plt.plot (gi, T4, color = 'black',  linewidth = 2, linestyle = 'dashed',  label = r"$c_{\beta\,k} = 0.1$")
plt.plot (gi, T2, color = 'black',  linewidth = 2, linestyle = 'dashdot', label = r"$c_{\beta\,k} = 0.5$")
plt.plot (gi, T1, color = 'black',  linewidth = 2, linestyle = 'solid',   label = r"$c_{\beta\,k} = 1.0$")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
#plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_k-Q_{E\,k}$', fontsize = "15")
plt.ylabel (r'$\hat{T}_k$',    fontsize = "15")
plt.legend (fontsize = "12");

plt.tight_layout ()

#plt.show () 
plt.savefig ("Fig8_4.pdf")
