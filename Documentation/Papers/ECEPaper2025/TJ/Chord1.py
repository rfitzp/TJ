# Chord.py

# Plots quantities along tilted central chord.
# User prompted for rational surface number and toroidal gridpoint.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ1.nc'
ds   = nc.Dataset(fn)
dTec = ds['Te_eq']
dTeO = ds['Te_1^O']
dTeX = ds['Te_2^X']
x    = ds['R_eq']
rres = ds['r_res']
Lres = np.asarray(ds['R_res'])
pn   = dTec.shape[2]

fn1   = 'TJ2.nc'
ds1   = nc.Dataset(fn1)
dTec1 = ds1['Te_eq']
dTeO1 = ds1['Te_1^O']
dTeX1 = ds1['Te_2^X']
x1    = ds1['R_eq']
rres1 = ds1['r_res']
Lres1 = np.asarray(ds1['R_res'])

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Quantities Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (2, 2, 2)

plt.xlim(1., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

k = 0
n1 = 0

for xx in Lres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (x, dTec[k,:,n1], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.plot    (x, dTeO[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')  
plt.plot    (x, dTeX[k,:,n1], color = 'red',   linewidth = 2,   linestyle = 'solid')  
plt.axhline (0.,              color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$',   fontsize = "17")
plt.ylabel (r'$T_e(eV)$',        fontsize = "17")

plt.subplot (2, 2, 1)

plt.xlim(1., x1[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

k = 0
n1 = 0

for xx in Lres1:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (x1, dTec1[k,:,n1], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.plot    (x1, dTeO1[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')  
plt.plot    (x1, dTeX1[k,:,n1], color = 'red',   linewidth = 2,   linestyle = 'solid')  
plt.axhline (0.,                color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$',   fontsize = "17")
plt.ylabel (r'$T_e(eV)$',        fontsize = "17")

plt.subplot (2, 2, 4)

plt.xlim(1., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

k = 0
n1 = 16

for xx in Lres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (x, dTec[k,:,n1], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.plot    (x, dTeO[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (x, dTeX[k,:,n1], color = 'red',   linewidth = 2,   linestyle = 'solid') 
plt.axhline (0.,              color = 'black', linewidth = 1.5, linestyle = 'dotted')

ymin, ymax = plt.ylim()

plt.xlabel (r'$R_\omega/R_0$',   fontsize = "17")
plt.ylabel (r'$T_e(eV)$',        fontsize = "17")

plt.subplot (2, 2, 3)

plt.xlim(1., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

k = 0
n1 = 8

for xx in Lres1:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (x1, dTec1[k,:,n1], color = 'green', linewidth = 2,   linestyle = 'solid')
plt.plot    (x1, dTeO1[k,:,n1], color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.plot    (x1, dTeX1[k,:,n1], color = 'red',   linewidth = 2,   linestyle = 'solid')  
plt.axhline (0.,                color = 'black', linewidth = 1.5, linestyle = 'dotted')

ymin, ymax = plt.ylim()

plt.xlabel (r'$R_\omega/R_0$',   fontsize = "17")
plt.ylabel (r'$T_e(eV)$',        fontsize = "17")

plt.tight_layout ()

#plt.show ()    
plt.savefig("Fig17.pdf")
