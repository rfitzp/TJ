# Straight.py

# Script plots Stage2 straight coordinate system

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = '../../Outputs/Flux/Stage1.nc'
ds = nc.Dataset(fn)

para = ds['Parameters']
rl   = ds['RLIM']
zl   = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi  = para[5]
rax  = para[6]
zax  = para[7]

aspect = (zhi - zlow) / (rrgt - rlft)

fn1  = '../../Outputs/Flux/Stage2.nc'
ds1  = nc.Dataset(fn1)
rb   = ds1['RBPTS']
zb   = ds1['ZBPTS']
rrst = ds1['RRst']
zzst = ds1['ZZst']

RR = np.asarray(rrst)
ZZ = np.asarray(zzst)          

fig = plt.figure (figsize = (7.0 / aspect ** 0.5, 7.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title ("FLUX: r, theta Coordinate System")
plt.rc ('xtick', labelsize = 12) 
plt.rc ('ytick', labelsize = 12) 

plt.subplot (1, 1, 1)

nr  = RR.shape[0]
for i in range(0, nr, 6):
    rr1 = RR[i]
    zz1 = ZZ[i]
    plt.plot (rr1, zz1, color = 'blue', linewidth = 0.5)
nth = RR.shape[1]
for j in range(0, nth, 6):
    rr1 = RR[:, j]
    zz1 = ZZ[:, j]
    plt.plot (rr1, zz1, color = 'black', linewidth = 0.5)
rr1 = RR[:, 0]
zz1 = ZZ[:, 0]
plt.plot (rr1, zz1, color = 'red', linewidth = 1)
rr1 = RR[:, int((nth+1)/2)]
zz1 = ZZ[:, int((nth+1)/2)]
plt.plot (rr1, zz1, color = 'red', linewidth = 1)
rr1 = RR[:, int((nth+1)/4)]
zz1 = ZZ[:, int((nth+1)/4)]
plt.plot (rr1, zz1, color = 'red', linewidth = 1)
rr1 = RR[:, int(3*(nth+1)/4)]
zz1 = ZZ[:, int(3*(nth+1)/4)]
plt.plot (rr1, zz1, color = 'red', linewidth  = 1)   

plt.plot (rax, zax, 'kx',          markersize = 4)
plt.plot (rb, zb, color = 'green', linewidth  = 1)

plt.xlabel ('$R/R_0$', fontsize = 12)
plt.ylabel ("$Z/R_0$", fontsize = 12, labelpad=-9)

plt.show()
