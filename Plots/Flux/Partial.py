# Partial.py

# Script plots Stage2 partial derivatives

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = '../../Outputs/Flux/Stage2.nc'
ds = nc.Dataset(fn)

r     = ds['r']
psi   = np.asarray(ds['PsiN'])
theta = ds['theta']
dRdr  = ds['dRdr']
dRdth = ds['dRdth']
dZdr  = ds['dZdr']
dZdth = ds['dZdth']
Jac   = ds['Jac']
Jax   = ds['Jax']

fig = plt.figure (figsize = (12.0, 7.0))
fig.canvas.manager.set_window_title ("FLUX: Partial Derivatives")

XX, YY = np.meshgrid (psi, theta, indexing = 'ij')

RR    = np.asarray (dRdr)
rrmin = np.min (RR)
rrmax = np.max (RR)
levrr = np.linspace (rrmin, rrmax, 80)

RT    = np.asarray (dRdth)
rtmin = np.min (RT)
rtmax = np.max (RT)
levrt = np.linspace (rtmin, rtmax, 40)

ZR    = np.asarray (dZdr)
zrmin = np.min (ZR)
zrmax = np.max (ZR)
levzr = np.linspace (zrmin, zrmax, 80)

ZT    = np.asarray (dZdth)
ztmin = np.min (ZT)
ztmax = np.max (ZT)
levzt = np.linspace (ztmin, ztmax, 40)

JT    = np.asarray (Jac)
jtmin = np.min (JT)
jtmax = np.max (JT)
levjt = np.linspace (jtmin, jtmax, 40)
JX    = np.asarray (Jax)

plt.subplot (3,  2, 1)

plt.contour (XX, YY, RR, levrr, linewidths = [1.0])

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$\\partial R/\\partial r$")

plt.subplot (3, 2, 2)

plt.contour (XX, YY, RT, levrt, linewidths = [1.0])

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$(\\partial R/\\partial\\theta)/r$")

plt.subplot (3, 2, 3)

plt.contour (XX, YY, ZR, levrr, linewidths=[1.0])

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$\\partial Z/\\partial r$")

plt.subplot (3, 2, 4)

plt.contour (XX, YY, ZT, levrt, linewidths=[1.0])

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$(\\partial Z/\\partial\\theta)/r$")

plt.subplot (3, 2, 5)

plt.contour (XX, YY, JT, levjt, linewidths=[1.0]) 

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$J$")

plt.subplot (3, 2, 6)

plt.contour (XX, YY, JX, levjt, linewidths=[1.0])

plt.xlabel ('$\\Psi_N$')
plt.ylabel ("$\\theta/\\pi$")
plt.title ("$J$ (analytic)")

plt.tight_layout(pad=0.5)

plt.show()
