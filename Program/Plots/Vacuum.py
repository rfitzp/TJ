# Vacuum.py

# Visualizes vacuum matrices, P_m^m' and R_m^m'

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn    = 'TJ.nc'
ds    = nc.Dataset(fn)
pvacr = ds['Pvac_r']
pvaci = ds['Pvac_i']
rvacr = ds['Rvac_r']
rvaci = ds['Rvac_i']

fig = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize=10) 
plt.rc ('ytick', labelsize=10)

plt.subplot (2, 2, 1)

prrax = plt.matshow (pvacr, fignum=0)
plt.colorbar (prrax)
plt.title (r"$Re(P_m^m)$")

plt.subplot (2, 2, 2)

priax = plt.matshow (pvaci, fignum=0)
plt.colorbar (priax)
plt.title (r"$Im(P_m^m)$")

plt.subplot (2, 2, 3)

rrrax = plt.matshow (rvacr, fignum=0)
plt.colorbar (rrrax)
plt.title (r"$Re(R_m^m)$")

plt.subplot (2, 2, 4)

rriax = plt.matshow (rvaci, fignum=0)
plt.colorbar (rriax)
plt.title (r"$Im(R_m^m)$")

plt.tight_layout ()

plt.show ()
