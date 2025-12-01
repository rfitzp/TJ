# RIii regime

# Qe    = 0.15
# Qi    = -0.15
# cbeta = 1.
# D     = 1.
# P     = 1.e-4
# tau   = -Qe/Qi = 1.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

fn = 'FourField.nc'
ds  = nc.Dataset(fn)
gi  = - np.asarray(ds['g_i'])
d3r = np.asarray(ds['Delta3_r'])
d3i = - np.asarray(np.asarray(ds['Delta3_i']))
d4r = np.asarray(ds['Delta4_r'])
d4i = - np.asarray(np.asarray(ds['Delta4_i']))

df =  pd.read_csv("RIii_Vz.txt", delim_whitespace=True, skiprows=1)

Q   = - np.asarray(df.iloc[:,0])
D4r = df.iloc[:,2]
D4i = - np.asarray(df.iloc[:,3])

fontsize = 17
    
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (2, 1, 1)

plt.xlim (gi[0], gi[-1])

plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.15,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.15, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, d3r, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Re($\hat{\Delta}_3$)")
plt.plot (gi, d4r, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Re($\hat{\Delta}_4$)")
plt.plot (Q,  D4r, color = 'cyan',  linewidth = 2, linestyle = 'dotted', label = r"SLAYER")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$', fontsize = fontsize)
plt.legend (fontsize = fontsize);

plt.subplot (2, 1, 2)

plt.xlim (gi[0], gi[-1])

plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.15,  color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (-0.15, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (gi, d3i, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = r"Im($\hat{\Delta}_3$)")
plt.plot (gi, d4i, color = 'red',   linewidth = 2, linestyle = 'solid',  label = r"Im($\hat{\Delta}_4$)")
plt.plot (Q,  D4i, color = 'cyan',  linewidth = 2, linestyle = 'dotted', label = r"SLAYER")

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$Q_E$', fontsize = fontsize)
plt.legend (fontsize = fontsize);

plt.tight_layout ()

#plt.show () 
plt.savefig ("Figure3.pdf")
