import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

def f(x, d, gammaw):

    if (d == 0):
        return x - gammaw
    else:
        return math.sqrt(x/d) * math.tanh (math.sqrt(x*d)) - gammaw

dd = np.linspace (0., 1., 100)

gg1 = []
gg2 = []
gg3 = []
gg4 = []

for d in dd:

    gammaw = 0.5
    g = brentq(lambda x: f(x, d, gammaw), 0, 2*gammaw)
    gg1.append (g)

    gammaw = 1.0
    g = brentq(lambda x: f(x, d, gammaw), 0, 2*gammaw)
    gg2.append (g)

    gammaw = 1.5
    g = brentq(lambda x: f(x, d, gammaw), 0, 2*gammaw)
    gg3.append (g)

    gammaw = 2.0
    g = brentq(lambda x: f(x, d, gammaw), 0, 4*gammaw)
    gg4.append (g)

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Vertical Code: Resistive Wall Mode Growth-Rate')
plt.rc('xtick', labelsize = 20) 
plt.rc('ytick', labelsize = 20)

plt.subplot(1, 1, 1)

plt.xlim (0., 1.)

plt.plot    (dd, gg1, color = 'red',   linewidth = 2,   linestyle = 'solid', label = "$\hat{\gamma}_{thin} = 0.5$")
plt.plot    (dd, gg2, color = 'green', linewidth = 2,   linestyle = 'solid', label = "$\hat{\gamma}_{thin} = 1.0$")
plt.plot    (dd, gg3, color = 'blue',  linewidth = 2,   linestyle = 'solid', label = "$\hat{\gamma}_{thin} = 1.5$")
plt.plot    (dd, gg4, color = 'cyan',  linewidth = 2,   linestyle = 'solid', label = "$\hat{\gamma}_{thin} = 2.0$")
#plt.axhline (0.,     color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\delta_w$',     fontsize = "20")
plt.ylabel (r"$\hat{\gamma}$", fontsize = "20")
plt.legend (fontsize = "20")

plt.tight_layout ()

#plt.show ()
plt.savefig ("Fig10.pdf")
    
