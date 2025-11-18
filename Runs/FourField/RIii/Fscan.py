# RIii regime

# Qe    = -0.15
# Qi    = +0.15
# cbeta = 1.
# D     = 1.
# P     = 1.e-4
# tau   = -Qe/Qi = 1.

# Delta = 2.124 * [i (Q - Qi)]^0.25 [i (Q - Qe)]

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

def F_r(Q):

    Qe     = 0.15
    Qi     = -0.15
    cbeta  = 1.
    D      = 1.
    P      = 1.e-3
    tau    = 1.

    i = complex (0, 1)
    
    return 2.124 * (pow(i * (Q - Qi), 0.25) * i * (Q - Qe)).real

def F_i(Q):

    Qe     = 0.15
    Qi     = -0.15
    cbeta  = 1.
    D      = 1.
    P      = 1.e-3
    tau    = 1.

    i = complex (0, 1)
    
    return 2.124 * (pow(i * (Q - Qi), 0.25) * i * (Q - Qe)).imag

fn = 'FourField.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = ds['g_i']
d3r = np.asarray(ds['Delta3_r'])
d3i = np.asarray(ds['Delta3_i'])
d4r = np.asarray(ds['Delta4_r'])
d4i = np.asarray(ds['Delta4_i'])

D_r = []
D_i = []
D_m = []
for Q in gi:
    D_r.append(F_r(Q))
    D_i.append(F_i(Q))
    D_m.append((F_r(Q)*F_r(Q) + F_i(Q)*F_i(Q))**0.5)

df =  pd.read_csv("RIii_Vz.txt", delim_whitespace=True, skiprows=1)

Q   = df.iloc[:,0]
D4r = df.iloc[:,2]
D4i = df.iloc[:,3]
    
fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, d3r, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Re($\hat{\Delta}_3$)")
plt.plot (gi, d4r, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Re($\hat{\Delta}_4$)")
plt.plot (Q,  D4r, color = 'cyan',  linewidth = 2, linestyle = 'dotted', label = r"slayer")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, d3i, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Im($\hat{\Delta}_3$)")
plt.plot (gi, d4i, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Im($\hat{\Delta}_4$)")
plt.plot (Q,  D4i, color = 'cyan',  linewidth = 2, linestyle = 'dotted', label = r"slayer")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
