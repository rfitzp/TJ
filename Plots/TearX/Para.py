# Para.py

# Plots layer parameters versus radius

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/TearX/TearX.nc'
ds = nc.Dataset(fn)
r  = ds['r']
QE = ds['Q_E']
Qe = ds['Q_e']
Qi = ds['Q_i']
D  = ds['D']
PE = ds['P_E']
Pp = ds['P_phi']
cb = ds['c_beta']
S  = ds['Scale']

rres = ds['rres']

x = input("r_max: ")
r_max = float (x)

r1  = []
QE1 = []
Qe1 = []
Qi1 = []
D1  = []
PE1 = []
Pp1 = []
cb1 = []
S1  = []

for rr, qqE, qqe, qqi, dd, ppe, ppp, ccb, ss in zip (r, QE, Qe, Qi, D, PE, Pp, cb, S):

    if (rr > r_max):
        r1.append(rr)
        QE1.append(qqE)
        Qe1.append(qqe)
        Qi1.append(qqi)
        D1.append(dd)
        PE1.append(ppe)
        Pp1.append(ppp)
        cb1.append(ccb)
        S1.append(ss)

Qep = [x + y for x, y in zip (QE1, Qe1)]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: Layer Parameters')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (r_max, 1.)

plt.plot    (r1, QE1, color = 'green',  linewidth = 2,   linestyle = 'solid', label = r"$Q_E$")
plt.plot    (r1, Qe1, color = 'red',    linewidth = 2,   linestyle = 'solid', label = r"$Q_e$")
plt.plot    (r1, Qi1, color = 'blue',   linewidth = 2,   linestyle = 'solid', label = r"$Q_i$")
plt.plot    (r1, Qep, color = 'cyan',   linewidth = 2,   linestyle = 'solid', label = r"$Q_{e\,\perp}$")
plt.axhline (0.,      color = 'black',  linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (r_max, 1.)
 
plt.plot    (r1, D1, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,   color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.ylabel (r"$D$", fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (r_max, 1.)
 
plt.plot    (r1, PE1, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$P_E$")
plt.plot    (r1, Pp1, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$P_\phi$")
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (r_max, 1.)
 
plt.plot    (r1, cb1,  color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',       fontsize = "15")
plt.ylabel (r'$c_\beta$', fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (r_max, 1.)
 
plt.plot    (r1, S1,  color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1.5, linestyle = 'dotted')
for rs in rres:
    plt.axvline (rs, color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$r$',   fontsize = "15")
plt.ylabel (r'Scale', fontsize = "15")

plt.tight_layout ()

plt.show ()    
