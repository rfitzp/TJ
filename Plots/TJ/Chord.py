# Chord.py

# Plots quantities along tilted central chord.
# User prompted for rational surface number and toroidal gridpoint.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/TJ/TJ.nc'
ds   = nc.Dataset(fn)
Tec  = ds['Te_eq']
Ted  = ds['Te_ece']
dTec = ds['dTe_eq']
dTed = ds['dTe_ece']
x    = ds['L_eq']
rres = ds['r_res']
Lres = np.asarray(ds['L_res'])
pn   = Tec.shape[2] 

nres = len(rres)
m    = input ("rational surface number (%d .. %d) ? " % (1, nres))
k    = int(m) - 1
nn   = input ("toroidal gridpoint (%d .. %d) ? " % (1, pn))
n1   = int(nn) - 1

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Quantities Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 2, 1)

#plt.xlim(x[0], x[-1])
plt.xlim(0., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

for xx in Lres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.plot    (x, dTec[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,              color = 'black', linewidth = 1.5, linestyle = 'dotted')

ymin, ymax = plt.ylim()

plt.xlabel (r'$x/R_0$',          fontsize = "15")
plt.ylabel (r'$\delta T_e(eV)$', fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim(0., x[-1])
#plt.xlim(x[0], x[-1])
plt.ylim (ymin, ymax)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

for xx in Lres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.plot    (x, dTed[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')    
plt.axhline (0.,              color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/R_0$',                   fontsize = "15")
plt.ylabel (r'$\delta T_{e\,\,ece}(eV)$', fontsize = "15")

plt.subplot (2, 2, 3)

#plt.xlim(x[0], x[-1])
plt.xlim(0., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

for xx in Lres:
    plt.axvline (xx, color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot    (x, Tec[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,             color = 'black', linewidth = 1.5, linestyle = 'dotted')

ymin, ymax = plt.ylim()

plt.xlabel (r'$x/R_0$',   fontsize = "15")
plt.ylabel (r'$T_e(eV)$', fontsize = "15")

plt.subplot (2, 2, 4)

#plt.xlim(x[0], x[-1])
plt.xlim(0., x[-1])
plt.ylim (ymin, ymax)

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

for xx in Lres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.plot    (x, Ted[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,             color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/R_0$',            fontsize = "15")
plt.ylabel (r'$T_{e\,\,ece}(eV)$', fontsize = "15")

plt.tight_layout ()

plt.show ()    

