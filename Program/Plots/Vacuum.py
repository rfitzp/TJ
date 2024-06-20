# Vacuum.py

# Visualizes vacuum solution matrices, P_m^m', Q_m^m', R_m^m', and S_m^m'.

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
pvacr = ds['Pvac_r']
pvaci = ds['Pvac_i']
qvacr = ds['Qvac_r']
qvaci = ds['Qvac_i']
rvacr = ds['Rvac_r']
rvaci = ds['Rvac_i']
svacr = ds['Svac_r']
svaci = ds['Svac_i']

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize=10) 
plt.rc ('ytick', labelsize=10)

plt.subplot (2, 4, 1)

prrax = plt.matshow (pvacr, fignum=0)
plt.colorbar (prrax)
plt.title (r"$Re(P_m^m)$")

plt.subplot (2, 4, 5)

priax = plt.matshow (pvaci, fignum=0)
plt.colorbar (priax)
plt.title (r"$Im(P_m^m)$")

plt.subplot (2, 4, 2)

qrrax = plt.matshow (qvacr, fignum=0)
plt.colorbar (qrrax)
plt.title (r"$Re(Q_m^m)$")

plt.subplot (2, 4, 6)

qriax = plt.matshow (qvaci, fignum=0)
plt.colorbar (qriax)
plt.title (r"$Im(Q_m^m)$")

plt.subplot (2, 4, 3)

rrrax = plt.matshow (rvacr, fignum=0)
plt.colorbar (rrrax)
plt.title (r"$Re(R_m^m)$")

plt.subplot (2, 4, 7)

rriax = plt.matshow (rvaci, fignum=0)
plt.colorbar (rriax)
plt.title (r"$Im(R_m^m)$")

plt.subplot (2, 4, 4)

srrax = plt.matshow (svacr, fignum=0)
plt.colorbar (srrax)
plt.title(r"$Re(S_m^m)$")

plt.subplot (2, 4, 8)

sriax = plt.matshow (svaci, fignum=0)
plt.colorbar (sriax)
plt.title (r"$Im(S_m^m)$")

plt.tight_layout ()

plt.show ()
