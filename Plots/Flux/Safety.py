# Safety.py

# Script plots Stage2 safety-factor profile

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = '../../Outputs/Flux/Stage2.nc'
ds = nc.Dataset (fn)

q  = ds['q']
qx = ds['q_x']
pn = ds['PsiN']

fig = plt.figure (figsize = (12.0, 6.))
fig.canvas.manager.set_window_title ("FLUX: Safety-Factor Profile")
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.xlim (0.0, 1.0)

plt.plot (pn, qx, 'b--', linewidth = 1.5)
plt.plot (pn, q,  'r--', linewidth = 1.)

plt.xlabel('$\\Psi_N$', fontsize = '15')
plt.ylabel('$q$',       fontsize = '15')

plt.show()
