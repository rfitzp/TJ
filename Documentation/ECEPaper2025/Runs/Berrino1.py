# Berinno.py

# Berrino algorithm for location of magnetic island via ece diagnostic
# User prompted for rational surface number.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn   = 'TJ1.nc'
ds   = nc.Dataset(fn)
Ted  = np.asarray(ds['Te_ece'])
Rres = np.asarray(ds['R_res'])
th   = np.asarray(ds['itheta_eq'])
R    = np.asarray(ds['R_eq'])
ix   = Ted.shape[1] 
pn   = Ted.shape[2]

RR = R - R/th

fn3   = 'TJ3.nc'
ds3   = nc.Dataset(fn3)
Ted3  = np.asarray(ds3['Te_ece'])

fn5   = 'TJ5.nc'
ds5   = nc.Dataset(fn5)
Ted5  = np.asarray(ds5['Te_ece'])

fn2   = 'TJ2.nc'
ds2   = nc.Dataset(fn2)
Ted2  = np.asarray(ds2['Te_ece'])
Rres2 = np.asarray(ds2['R_res'])

fn4   = 'TJ4.nc'
ds4   = nc.Dataset(fn4)
Ted4  = np.asarray(ds4['Te_ece'])

fn6   = 'TJ6.nc'
ds6   = nc.Dataset(fn6)
Ted6  = np.asarray(ds6['Te_ece'])

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
            sum = sum  + (Ted[k, i - off, p] -  Ted[k, i + off, p])**2 /4./off/off

    ber.append(sum)

ber3 = []
for i in range (ignore):
    ber3.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (Ted3[k, i - off, p] -  Ted3[k, i + off, p])**2 /4./off/off

    ber3.append(sum)

ber5 = []
for i in range (ignore):
    ber5.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (Ted5[k, i - off, p] -  Ted5[k, i + off, p])**2 /4./off/off

    ber5.append(sum)

ber2 = []
for i in range (ignore):
    ber2.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (Ted2[k, i - off, p] -  Ted2[k, i + off, p])**2 /4./off/off

    ber2.append(sum)            

ber4 = []
for i in range (ignore):
    ber4.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (Ted4[k, i - off, p] -  Ted4[k, i + off, p])**2 /4./off/off

    ber4.append(sum)            

ber6 = []
for i in range (ignore):
    ber6.append(0.)

for i in range (ignore, ix):
    sum = 0.

    if (i - off > 0 and i + off < ix-1):
        for p in range (pn):
            sum = sum  + (Ted6[k, i - off, p] -  Ted6[k, i + off, p])**2 /6./off/off

    ber6.append(1.5*sum)            
    
fig = plt.figure (figsize = (12.0, 6.0))
fig.canvas.manager.set_window_title (r'TJ Code: Berrino Algorithm')
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (1, 2, 2)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (RR, ber,  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (RR, ber3, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (RR, ber5, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega\,(1-\theta_w)/R_0$', fontsize = "17")
plt.legend (fontsize = "15")

plt.subplot (1, 2, 1)

plt.xlim(1., R[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))
    
for xx in Rres2:
    plt.axvline (xx, color = 'black',  linewidth = 1.5, linestyle = 'dashed')

plt.plot    (RR, ber2, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$W/a=0.10$")
plt.plot    (RR, ber4, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$W/a=0.05$")
plt.plot    (RR, ber6, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$W/a=0.01$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$R_\omega\,(1-\theta_w)/R_0$', fontsize = "17")
plt.legend (fontsize = "15")
                        
plt.tight_layout ()

#plt.show ()    
plt.savefig ("Berrino1.pdf")
