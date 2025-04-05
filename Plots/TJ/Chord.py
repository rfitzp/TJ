# Chord.py

# Plots perturbed quantities along tilted central chord

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn1 = '../../Outputs/Equilibrium/Equilibrium.nc'
ds1 = nc.Dataset(fn1)
BR  = np.asarray(ds1['BR_eq'])
ne  = np.asarray(ds1['ne_eq'])

fn   = '../../Outputs/TJ/TJ.nc'
ds   = nc.Dataset(fn)
bRc  = np.asarray(ds['b_R_eq_cos'])
bRs  = np.asarray(ds['b_R_eq_sin'])
dnec = np.asarray(ds['dne_eq_cos'])
dnes = np.asarray(ds['dne_eq_sin'])
dTec = ds['dTe_eq_cos']
dTes = ds['dTe_eq_sin']
x    = ds['L_eq']
rres = ds['r_res']

nres = len(rres)
print ("rational surface = (%d .. %d)" % (1, nres))
m = input ("rational surface number ? ")
k = int(m) - 1

dbnc = bRc[k,:] * ne + BR * dnec[k,:]
dbns = bRs[k,:] * ne + BR * dnes[k,:]

intc = np.trapz (dbnc, x)
ints = np.trapz (dbns, x)

print ("cos = %10.3e sin = %10.3e" % (intc, ints))

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Perturbed Quantities Along Tilted Central Chord')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (3, 1, 1)

plt.xlim(0., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (x, dTec[k,:], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (x, dTes[k,:], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,                  color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/R_0$',           fontsize = "15")
plt.ylabel (r'$\delta T_e(eV)$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim(0., x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (x, dnec[k,:], color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (x, dnes[k,:], color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (0.,           color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/R_0$',               fontsize = "15")
plt.ylabel (r'$\delta n_e(m^{-3})$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim(x[0], x[-1])

plt.gca().ticklabel_format (axis='y', style='sci', scilimits=(0, 0))

plt.plot    (x, dbnc, color = 'blue',   linewidth = 2,   linestyle = 'solid',  label = 'cos')
plt.plot    (x, dbns, color = 'green',  linewidth = 2,   linestyle = 'solid',  label = 'sin')
plt.axhline (intc,    color = 'blue',   linewidth = 1.5, linestyle = 'dotted')
plt.axhline (ints,    color = 'green',  linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.,     color = 'black',   linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$x/R_0$',                                 fontsize = "15")
plt.ylabel (r'$\delta (B_\parallel\,n_e)(T\,m^{-3})$',  fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    

