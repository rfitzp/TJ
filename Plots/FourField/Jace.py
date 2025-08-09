# Fscan.py

# Plots real and imaginary parts of Delta3, Delta4, DeltaJ1, and DeltaJ2 versus imagninary part of g

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
j1r = ds['DeltaJ1_r']
j1i = ds['DeltaJ1_i']
j2r = np.asarray(ds['DeltaJ2_r'])
j2i = - np.asarray(ds['DeltaJ2_i'])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Delta3 and Delta4')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 1, 1)

#plt.xlim (gi[0], gi[-1])
plt.xlim (-0.15, -0.05)
plt.ylim (-0.1,   0.1)

plt.plot (gi, d3r, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_3$)")
plt.plot (gi, d4r, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_4$)")
plt.plot (gi, j1r, color = 'green', linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_{J1}$)")
plt.plot (gi, j2r, color = 'cyan',  linewidth = 2, linestyle = 'solid', label = r"Re($\Delta_{J2}$)")

plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Im(g)$', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (2, 1, 2)

#plt.xlim (gi[0], gi[-1])
plt.xlim (-0.15, -0.05)
plt.ylim (-0.1,   0.1)

plt.plot (gi, d3i, color = 'blue', linewidth = 2, linestyle = 'solid',  label = r"Im($\Delta_3$)")
plt.plot (gi, d4i, color = 'red',  linewidth = 2, linestyle = 'solid',  label = r"Im($\Delta_4$)")
plt.plot (gi, j1i, color = 'green', linewidth = 2, linestyle = 'solid', label = r"Im($\Delta_{J1}$)")
plt.plot (gi, j2i, color = 'cyan',  linewidth = 2, linestyle = 'solid', label = r"Im($\Delta_{J2}$)")


plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[2], color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-In[3], color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Im(g)$', fontsize = "15")
plt.legend (fontsize = "15");

plt.tight_layout ()

plt.show () 
#plt.savefig ("Fscan.pdf")
