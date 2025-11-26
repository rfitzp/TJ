# cbeta.py

# Plots Q_E values at which Im(Delta) = 0 versus cbeta

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df1 = pd.read_csv("HRi.txt", delim_whitespace=True, skiprows=1)

c1 = df1.iloc[:,0] 
q1 = df1.iloc[:,1] - 1.

df2 = pd.read_csv("HRii.txt", delim_whitespace=True, skiprows=1)

c2 = df2.iloc[:,0] 
q2 = df2.iloc[:,1] - 1.

df3 = pd.read_csv("VRii.txt", delim_whitespace=True, skiprows=1)

c3 = df3.iloc[:,0] 
q3 = df3.iloc[:,1] - 0.15

fontsize = 17

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (1, 1, 1)

plt.xlim (0., 1.0)

plt.plot (c1, q1, color = 'red',   linewidth = 2, linestyle = 'solid', label = "HRi")
plt.plot (c2, q2, color = 'green', linewidth = 2, linestyle = 'solid', label = "HRii")
plt.plot (c3, q3, color = 'blue',  linewidth = 2, linestyle = 'solid', label = "VRi")

plt.plot ([0.,c1[0]], [0., q1[0]], color = 'red',   linewidth = 2, linestyle = 'solid')
plt.plot ([0.,c2[0]], [0., q2[0]], color = 'green', linewidth = 2, linestyle = 'solid')
plt.plot ([0.,c3[0]], [0., q3[0]], color = 'blue',  linewidth = 2, linestyle = 'solid')

plt.axhline (0., color = 'black',   linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$c_\beta$',              fontsize = fontsize)
plt.ylabel (r'$(Q_{E} + Q_e)_{crit}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

plt.show () 
#plt.savefig ("Fscan.pdf")
