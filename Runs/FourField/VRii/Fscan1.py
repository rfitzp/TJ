# VRii regime

# Qe    = -0.15
# Qi    = +0.15
# cbeta = 1.
# D     = 0.1
# P     = 1.
# tau   = -Qe/Qi = 1.

# Delta = 2.104 * [i (Q - Qe)] P^1/6

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

def F_r(Q):

    Qe     = 0.15
    Qi     = -0.15
    cbeta  = 1.
    D      = 0.1
    P      = 1.
    tau    = 1.

    return 0.
 
def F_i(Q):

    Qe     = 0.15
    Qi     = -0.15
    cbeta  = 1.
    D      = 0.1
    P      = 1.
    tau    = 1.
    
    return 2.104 * (Q - Qe) * pow (P, 1./6.)

df  = pd.read_csv("slayer_Qscan_n1.out", skiprows=1, delim_whitespace=True, header=None)

fn = 'FourField.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = ds['g_i']
d3r = np.asarray(ds['Delta3_r'])
d3i = np.asarray(ds['Delta3_i'])
d4r = np.asarray(ds['Delta4_r'])
d4i = np.asarray(ds['Delta4_i'])

d3m = (d3r*d3r + d3i*d3i)**0.5
d4m = (d4r*d4r + d4i*d4i)**0.5

D_r = []
D_i = []
D_m = []
for Q in gi:
    D_r.append(F_r(Q))
    D_i.append(F_i(Q))
    D_m.append((F_r(Q)*F_r(Q) + F_i(Q)*F_i(Q))**0.5)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.plot (gi, d3r, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Re($\Delta_3$)")
plt.plot (gi, d4r, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Re($\Delta_4$)")
plt.plot (gi, D_r, color = 'green', linewidth = 2, linestyle = 'dotted', label = r"Analytic")

plt.plot (df[0].values, df[2].values, color = 'cyan', linewidth = 2, linestyle = 'solid', label = "slayer")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])
plt.gca().invert_xaxis()

plt.plot (gi, d3i, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Im($\Delta_3$)")
plt.plot (gi, d4i, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Im($\Delta_4$)")
plt.plot (gi, D_i, color = 'green', linewidth = 2, linestyle = 'dotted', label = r"Analytic")

plt.plot (df[0].values, df[3].values, color = 'cyan', linewidth = 2, linestyle = 'solid', label = "slayer")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
#plt.savefig ("VRii.pdf")

with open("VRii.txt", "w") as f:
    f.write ("RIii: Qe = 0.15, Qi = -0.15, cbeta = 1, D = 1, P = 1.e-3: Columns are Q, Delta3_r, Delta3_i, Delta4_r, Delta4_i, Delta_analytic_r, Delta_analytic_i\n")
    for gi, x3, y3, x4, y4, xa, ya in zip (gi, d3r, d3i, d4r, d4i, D_r, D_i):
        f.write ("%10.4e  %10.3e %10.3e  %10.3e %10.3e  %10.3e %10.3e\n" % (gi, x3, y3, x4, y4, xa, ya))

f.close ()        
