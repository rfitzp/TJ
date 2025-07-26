# Berinno1.py

# Modified Berrino algorithm for location of magnetic island via ece diagnostic
# User prompted for rational surface number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = '../../Outputs/TJ/TJ.nc'
ds   = nc.Dataset(fn)
TeO  = np.asarray(ds['Te_1^O'])
TeX  = np.asarray(ds['Te_2^X'])
DO   = np.asarray(ds['Delta_1^O'])
DX   = np.asarray(ds['Delta_2^X'])
Rres = np.asarray(ds['R_res'])
Ores = np.asarray(ds['O_res'])
Xres = np.asarray(ds['X_res'])
R    = np.asarray(ds['R_eq'])
ix   = TeO.shape[1] 
pn   = TeO.shape[2]

nres = len(Rres)
m    = input ("rational surface number (%d .. %d) ? " % (1, nres))
k    = int(m) - 1
mm   = input ("offset ? ")
off  = int(mm)

RO = R - DO
RX = R - DX

RO -= (Rres[k] - Ores[k])
RX -= (Rres[k] - Ores[k])

ignore = 1000

berO = []
for i in range (ignore):
    berO.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO[k, i - off, p] -  TeO[k, i + off, p])**2 /4./off/off

    berO.append(sum)

berX = []
for i in range (ignore):
    berX.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX[k, i - off, p] -  TeX[k, i + off, p])**2 /4./off/off

    berX.append(sum)        

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Berrino Algorithm')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

for xx in Ores:
    plt.axvline (xx, color = 'red',  linewidth = 1.5, linestyle = 'dotted')

for xx in Xres:
    plt.axvline (xx, color = 'blue', linewidth = 1.5, linestyle = 'dotted')
    
for xx in Rres:
    plt.axvline (xx, color = 'black', linewidth = 1.5, linestyle = 'dashed')

plt.plot    (RO, berO, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = 'O-mode')
plt.plot    (RX, berX, color = 'red',   linewidth = 2,   linestyle = 'solid', label = 'X-mode')
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$(R-\Delta - \delta\,W)/R_0$', fontsize = "15")
plt.legend (fontsize = "15")
            
plt.tight_layout ()

plt.show ()    

