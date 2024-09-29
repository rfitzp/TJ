# Solutions.py

# Plots torque and shielding curves associated with kth rational surface in plasma
# User prompted for k.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Layer/layer.nc'
ds = nc.Dataset(fn)
w  = ds['omega_r']
dr = ds['Deltar']
di = ds['Deltai']
T  = ds['T_res']
X  = ds['Xi_res']

fn1 = '../../Outputs/TJ/Tj.nc'
ds1 = nc.Dataset(fn1)
Qe  = ds1['Qe']
QE  = ds1['QE']
Qi  = ds1['Qi']
tau = ds1['tau']          

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Torques and Shielding Factors')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

nres = T.shape[0]
print ("rational surfaces = (%d .. %d)" % (1, nres))
m = input ("k ? ")
k = int(m) - 1

ww = np.asarray(w[k,:])/1.e3;

plt.subplot (3, 1, 1)

plt.xlim (ww[0], ww[-1])

plt.plot (ww, dr[k,:], color = 'blue', linewidth = 2, linestyle = 'solid', label = r"Re($\Delta$)")
plt.plot (ww, di[k,:], color = 'red',  linewidth = 2, linestyle = 'solid', label = r"Im($\Delta$)")

plt.axhline (0.,                        color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k]+QE[k])/tau[k]/1.e3, color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k])/tau[k]/1.e3,       color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qi[k]+QE[k])/tau[k]/1.e3, color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$\omega-\omega_E$ (kHz)', fontsize = "15")
plt.legend (fontsize = "15");

plt.subplot (3, 1, 2)

plt.xlim (ww[0], ww[-1])

plt.plot (ww, T[k,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0.,                        color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k]+QE[k])/tau[k]/1.e3, color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k])/tau[k]/1.e3,       color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qi[k]+QE[k])/tau[k]/1.e3, color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$\omega-\omega_E$ (kHz)', fontsize = "15")
plt.ylabel (r"$\delta T$",              fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (ww[0], ww[-1])

plt.plot (ww, X[k,:], color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0.,                        color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k]+QE[k])/tau[k]/1.e3, color = 'red',   linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qe[k])/tau[k]/1.e3,       color = 'green', linewidth = 1.5, linestyle = 'dotted')
plt.axvline ((Qi[k]+QE[k])/tau[k]/1.e3, color = 'blue',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\omega-\omega_E$ (kHz)', fontsize = "15")
plt.ylabel (r"$\Xi$",                   fontsize = "15")

plt.tight_layout ()

plt.show () 
