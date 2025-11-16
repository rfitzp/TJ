# Psinu.py

# Plots PSI versus q95 for nu scan

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("../../Outputs/TearX/Scannu.txt", delim_whitespace = True, header = None)

q95 = df.iloc[:,4]
dlt = df.iloc[:,8]
psi = df.iloc[:,9]

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: nu Scan')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.ylim (0.9, 1.001)

plt.scatter (q95, psi, c = 'blue', s = 10)

plt.axhline (1.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axhline (0.95,  color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$q_{95}$', fontsize = "15")
plt.ylabel (r"$\Psi$",   fontsize = "15")

plt.tight_layout ()

plt.show ()    
