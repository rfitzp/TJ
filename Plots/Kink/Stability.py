import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '../../Outputs/Kink/Kink.nc'
ds = nc.Dataset(fn)
qa = ds['q_xa']
nu = ds['nu']
p  = ds['para']

mpol = p[1]
ntor = p[2]

qs = mpol/ntor

fontsize = 15

fig = plt.figure (figsize = (10.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

numax = 1.05*max(np.asarray(nu))

#plt.xlim (1., 8.1)
plt.ylim (1., numax)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.plot (qa, nu, color = 'black',  linewidth = 2, linestyle = 'solid')

plt.xlabel (r'$q_a$', fontsize = fontsize)
plt.ylabel (r'$\nu$', fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
#plt.savefig("Figure9-17.pdf")
