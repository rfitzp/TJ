# Plot.py

# Plots aLL, aLS, aSL, and aSS versus nu

import math
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

def aLL(nu):

    g1 = gamma(1.  + nu/2.)
    g2 = gamma(1.5 + nu)
    g3 = gamma(0.5 + nu/2.)
    g4 = gamma(1.  + nu)

    return 2.**(1.+nu) * g1 * g2 /(1. + nu) /g3 /g4

def aSL(nu):

    g1 = gamma(1.  + nu/2.)
    g2 = gamma(1.5 + nu)
    g3 = gamma(0.5 + nu/2.)
    g4 = gamma(1.  + nu)

    return  - 2.**nu * nu * g3 * g2 /g1 /g4

def aLS(nu):

    g1 = gamma(0.5 - nu/2.)
    g2 = gamma(0.5 - nu)
    g3 = gamma(1.  - nu)
    g4 = gamma(1.  - nu/2.)

    return  - nu * g1 * g2 /2.**(1.+nu) /g3 /g4

def aSS(nu):

    g1 = gamma(0.5 - nu/2.)
    g2 = gamma(0.5 - nu)
    g3 = gamma(1.  - nu)
    g4 = gamma(1.  - nu/2.)

    return (1. + nu) * g4 * g2 /2.**nu /g1 /g3

xmin = -0.5
xmax = 0.5

x = np.linspace(xmin, xmax, 1000)

y1 = aLL(x)
y2 = aSL(x)
y3 = aLS(x)
y4 = aSS(x)

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.xlim (xmin, xmax)
plt.ylim (-5., 5.)
 
plt.plot    (x, y1, color = 'red',   linewidth = 2,   linestyle = 'solid', label = r'$a_{LL}$')
plt.plot    (x, y2, color = 'green', linewidth = 2,   linestyle = 'solid', label = r'$a_{SL}$')
plt.plot    (x, y3, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = r'$a_{LS}$')
plt.plot    (x, y4, color = 'cyan',  linewidth = 2,   linestyle = 'solid', label = r'$a_{SS}$')
plt.axhline (0.,    color = 'black', linewidth = 1.,  linestyle = 'dotted')
plt.axvline (0.,    color = 'black', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$\nu$',    fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    
