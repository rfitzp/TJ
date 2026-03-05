import matplotlib.pyplot as plt
import math
import netCDF4 as nc
import numpy as np
from scipy import integrate

def f (x, nu):
    return (1. - (1. - x)**(nu+1.))**2 /x

nu = []
li = []

for n in np.arange (1., 5., 1.e-3):
    
    res, err = integrate.quad (f, 1.e-3, 1., args = (n))

    nu.append (n)
    li.append (res)
 
fig = plt.figure (figsize=(8.0, 8.0))
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 

plt.subplot (1, 1, 1)

plt.grid ()
plt.xlim (1., 5.)
plt.ylim(0.9, 1.8)
plt.yticks([1, 1.2, 1.4, 1.6, 1.8])

plt.plot (nu, li, color = 'black',  linewidth = 2, linestyle = 'solid')

plt.xlabel(r'$\nu$', fontsize = "20")
plt.ylabel(r'$l_i$', fontsize = '20')

#plt.show()
plt.savefig("Figure9_19.pdf")
    

