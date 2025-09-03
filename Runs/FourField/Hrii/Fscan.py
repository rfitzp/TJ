# HRii regime

# Qe    = -1.
# Qi    = +1.
# cbeta = 1.
# D     = 10.
# P     = 0.1
# tau   = -Qe/Qi = 1.

# Delta = 2.124 * i (Q - Qe) cbeta^1/2 D^-1/2 

def F_r(Q):

    Qe     = 1.
    cbeta  = 1.
    D      = 10.
    P      = 0.1
    tau    = 1.

    return 0.

def F_i(Q):

    Qe     = 1.
    cbeta  = 1.
    D      = 10.
    P      = 0.1
    tau    = 1.

    return 2.124 * (Q - Qe) * cbeta**0.5 /D**0.5 

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'FourField.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = ds['g_i']
d3r = ds['Delta3_r']
d3i = ds['Delta3_i']
d4r = ds['Delta4_r']
d4i = ds['Delta4_i']

D_r = []
D_i = []
for Q in gi:
    D_r.append(F_r(Q))
    D_i.append(F_i(Q))

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

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
#plt.savefig ("HRii.pdf")

with open("Hrii.txt", "w") as f:
    f.write ("HRii: Qe = 1, Qi = -1, cbeta = 1, D = 10, P = 0.1: Columns are Q, Delta3_r, Delta3_i, Delta4_r, Delta4_i, Delta_analytic_r, Delta_analytic_i\n")
    for gi, x3, y3, x4, y4, xa, ya in zip (gi, d3r, d3i, d4r, d4i, D_r, D_i):
        f.write ("%10.4e  %10.3e %10.3e  %10.3e %10.3e  %10.3e %10.3e\n" % (gi, x3, y3, x4, y4, xa, ya))

f.close ()        
