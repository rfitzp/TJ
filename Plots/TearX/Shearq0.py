# Shearq0.py

# Plots magnetic shear versus q95 for q0 scan

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("../../Outputs/TearX/Scanq0.txt", delim_whitespace = True, header = None)

q95 = df.iloc[:,4]
ss  = df.iloc[:,7]  

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: q0 Scan')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.scatter (q95, ss, c = 'blue', s = 10)

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$q_{95}$', fontsize = "15")
plt.ylabel (r'$s$',      fontsize = "15")

plt.tight_layout ()

plt.show ()    
