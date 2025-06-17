# Fscan.py

# Plots real and imaginary parts of Delta3 and Delta4 versus imagninary part of g
# User prompted for k.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/FourField/FourField.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']
gi  = ds['g_i']
d3r = ds['Delta3_r']
d3i = ds['Delta3_i']
d4r = ds['Delta4_r']
d4i = ds['Delta4_i']

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])

plt.plot (gi, d3r, color = 'blue', linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_3$)")
plt.plot (gi, d4r, color = 'red',  linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_4$)")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Im(g)$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])

plt.plot (gi, d3i, color = 'blue', linewidth = 2, linestyle = 'solid', label = r"Im($\Delta_3$)")
plt.plot (gi, d4i, color = 'red',  linewidth = 2, linestyle = 'solid', label = r"Im($\Delta_4$)")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Im(g)$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
