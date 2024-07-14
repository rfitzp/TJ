# Shaping.py

# Plots shaping functions versus r
 
import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = 'Equilibrium.nc'
ds  = nc.Dataset(fn)
r   = ds['r']
g2  = ds['g_2']
Hn  = ds['Hn']
Hnp = ds['Hnp']
Vn  = ds['Vn']
Vnp = ds['Vnp']

Hnn = np.asarray(Hn)
ns   = Hnn.shape[0] - 1
ns   = 5
                      
fig = plt.figure(figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: Shaping Functions')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (ns, 2, 1)

plt.xlim (0., 1.)

plt.plot    (r, g2, color = 'blue',  linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$g_2$',     fontsize = "15")

plt.subplot (ns, 2, 2)
 
plt.xlim (0., 1.)

plt.plot    (r, Hn[1],  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$H_1$")
plt.plot    (r, Hnp[1], color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$H_1'$")
plt.axhline (0.,        color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize="15")
plt.ylabel (r'$H_1$',     fontsize = "15")

j = 3
for n in range (2, ns+1):

    plt.subplot (ns, 2, j)
    j = j + 1

    plt.xlim (0., 1.)

    plt.plot    (r, Hn[n],  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$H_%1d$"  % n)
    plt.plot    (r, Hnp[n], color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$H_%1d'$" % n)
    plt.axhline (0.,        color = 'black', linewidth = 1.5, linestyle = 'dotted')

    plt.xlabel (r'$\hat{r}$', fontsize="15")
    plt.ylabel (r'$H_{%2d}$'%(n), fontsize = "15")

    plt.subplot (ns, 2, j)
    j = j +1

    plt.xlim (0., 1.)

    plt.plot    (r, Vn[n],  color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r"$V_%1d$"  % n)
    plt.plot    (r, Vnp[n], color = 'red',   linewidth = 2,   linestyle = 'solid', label = r"$V_%1d'$" % n)
    plt.axhline (0.,        color = 'black', linewidth = 1.5, linestyle = 'dotted')

    plt.xlabel (r'$\hat{r}$', fontsize="15")
    plt.ylabel (r'$V_{%2d}$'%(n), fontsize = "15")

plt.tight_layout ()

plt.show ()    
