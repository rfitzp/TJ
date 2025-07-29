# Berinno.py

# Berrino algorithm for location of magnetic island via ece diagnostic
# User prompted for rational surface number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ1.nc'
ds   = nc.Dataset(fn)
TeO  = np.asarray(ds['Te_1^O'])
TeX  = np.asarray(ds['Te_2^X'])
Rres = np.asarray(ds['R_res'])
R    = ds['R_eq']
ix   = TeO.shape[1] 
pn   = TeO.shape[2]

fn3   = 'TJ3.nc'
ds3   = nc.Dataset(fn3)
TeO3  = np.asarray(ds3['Te_1^O'])
TeX3  = np.asarray(ds3['Te_2^X'])

fn5   = 'TJ5.nc'
ds5   = nc.Dataset(fn5)
TeO5  = np.asarray(ds5['Te_1^O'])
TeX5  = np.asarray(ds5['Te_2^X'])

fn2   = 'TJ2.nc'
ds2   = nc.Dataset(fn2)
TeO2  = np.asarray(ds2['Te_1^O'])
TeX2  = np.asarray(ds2['Te_2^X'])
Rres2 = np.asarray(ds2['R_res'])

fn4   = 'TJ4.nc'
ds4   = nc.Dataset(fn4)
TeO4  = np.asarray(ds4['Te_1^O'])
TeX4  = np.asarray(ds4['Te_2^X'])

fn6   = 'TJ6.nc'
ds6   = nc.Dataset(fn6)
TeO6  = np.asarray(ds6['Te_1^O'])
TeX6  = np.asarray(ds6['Te_2^X'])

k = 0
off = 1
ignore = 1000

ber = []
for i in range (ignore):
    ber.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO[k, i - off, p] -  TeO[k, i + off, p])**2 /4./off/off

    ber.append(sum)

ber3 = []
for i in range (ignore):
    ber3.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO3[k, i - off, p] -  TeO3[k, i + off, p])**2 /4./off/off

    ber3.append(sum)

ber5 = []
for i in range (ignore):
    ber5.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO5[k, i - off, p] -  TeO5[k, i + off, p])**2 /4./off/off

    ber5.append(sum)

ber2 = []
for i in range (ignore):
    ber2.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO2[k, i - off, p] -  TeO2[k, i + off, p])**2 /4./off/off

    ber2.append(sum)            

ber4 = []
for i in range (ignore):
    ber4.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO4[k, i - off, p] -  TeO4[k, i + off, p])**2 /4./off/off

    ber4.append(sum)            

ber6 = []
for i in range (ignore):
    ber6.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeO6[k, i - off, p] -  TeO6[k, i + off, p])**2 /6./off/off

    ber6.append(1.5*sum)            

berX = []
for i in range (ignore):
    berX.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX[k, i - off, p] -  TeX[k, i + off, p])**2 /4./off/off

    berX.append(sum)

berX3 = []
for i in range (ignore):
    berX3.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX3[k, i - off, p] -  TeX3[k, i + off, p])**2 /4./off/off

    berX3.append(sum)

berX5 = []
for i in range (ignore):
    berX5.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX5[k, i - off, p] -  TeX5[k, i + off, p])**2 /4./off/off

    berX5.append(sum)

berX2 = []
for i in range (ignore):
    berX2.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX2[k, i - off, p] -  TeX2[k, i + off, p])**2 /4./off/off

    berX2.append(sum)            

berX4 = []
for i in range (ignore):
    berX4.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX4[k, i - off, p] -  TeX4[k, i + off, p])**2 /4./off/off

    berX4.append(sum)            

berX6 = []
for i in range (ignore):
    berX6.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (TeX6[k, i - off, p] -  TeX6[k, i + off, p])**2 /6./off/off

    berX6.append(1.5*sum)            
    
fig = plt.figure (figsize = (10.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Berrino Algorithm')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (2, 2, 2)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, ber,  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (R, ber3, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (R, ber5, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$', fontsize = "17")
plt.ylabel (r'O-mode (a.u.)',   fontsize = "17")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 1)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres2:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, ber2, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (R, ber4, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (R, ber6, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$', fontsize = "17")
plt.ylabel (r'O-mode (a.u.)',   fontsize = "17")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, berX,  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (R, berX3, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (R, berX5, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$', fontsize = "17")
plt.ylabel (r'X-mode (a.u.)', fontsize = "17")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres2:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (R, berX2, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (R, berX4, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (R, berX6, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega/R_0$', fontsize = "17")
plt.ylabel (r'X-mode (a.u.)', fontsize = "17")
plt.legend (fontsize = "15")
                        
plt.tight_layout ()

#plt.show ()    
plt.savefig ("Fig18.pdf")
