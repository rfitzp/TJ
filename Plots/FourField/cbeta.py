# cbeta.py

# Plots Q_E values at which Im(Delta) = 0 versus cbeta

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

fn = '../../Outputs/FourField/FourField.nc'
ds  = nc.Dataset(fn)
In  = ds['InputParameters']

df = pd.read_csv("../../Outputs/FourField/cbetascan.txt", delim_whitespace=True, skiprows=1)

c = df.iloc[:,0] 
q = df.iloc[:,1] + In[2]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Q_E_crit versus c_beta')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

plt.xlim (0., 1.0)

plt.plot (c, q, color = 'blue', linewidth = 2, linestyle = 'solid')
#plt.plot ([0.,c[0]], [1.0,q[0]], color = 'blue', linewidth = 2, linestyle = 'solid')

#plt.axhline (1., color = 'black',   linewidth = 1.5, linestyle = 'dotted')
#plt.axhline (0.,     color = 'green', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$c_\beta$',          fontsize = "15")
plt.ylabel (r'$(Q_E)_{crit}$', fontsize = "15")

plt.tight_layout ()

plt.show () 
#plt.savefig ("Fscan.pdf")
