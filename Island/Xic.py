# Script to plot function that determines xic

import math
import numpy as np
import scipy.special as sp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

Nh = 8
Nb = 8

delta = 0.5

def csz (xi):

    sum = - delta*delta/2.
    
    for n in range (1, nh):
        sum += math.cos(n * xi) * (sp.Jn(n-1, n*delta*delta) - sp.Jn(n+1, n*delta*delta))/n
                
    return sum

def fun (xi, X):

    cs = math.cos (xi)

    return 1. - X*X - math.cos(xi) + 2.*math.sqrt (8.) * delta * X * cs - delta*delta * cs*cs

def xic (x):

    if (fun(0., x) * fun (math.pi, x) > 0.):
        return math.pi
    else:
        result = root_scalar(lambda xi: fun(xi, x), bracket=[0., math.pi], method='brentq')
        return result.root
    

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TJ Code: \xi_c(X)')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15)

plt.subplot (1, 1, 1)

#plt.xlim (x[0], x[-1])
#plt.ylim (0.,   1.)

xx = np.arange (-4., 4., 0.01)

f = []

for x in xx:
    f.append (xic (x) /math.pi)

plt.xlim (-4., 4.)

plt.plot    (xx, f,  color = 'blue', linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\xi/\pi$', fontsize = "15")
plt.ylabel (r'$F$',       fontsize = "15")

plt.tight_layout ()

plt.show () 
