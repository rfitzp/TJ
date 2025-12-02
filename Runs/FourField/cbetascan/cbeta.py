import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df1 = pd.read_csv("HRi.txt", delim_whitespace=True, skiprows=1)

c1 = df1.iloc[:,0] 
q1 = - df1.iloc[:,1] + 1.

df1a = pd.read_csv("HRia.txt", delim_whitespace=True, skiprows=1)

c1a = df1a.iloc[:,0] 
q1a = - df1a.iloc[:,1] + 1.

df2 = pd.read_csv("HRii.txt", delim_whitespace=True, skiprows=1)

c2 = df2.iloc[:,0] 
q2 = - df2.iloc[:,1] + 1.

df2a = pd.read_csv("HRiia.txt", delim_whitespace=True, skiprows=1)

c2a = df2a.iloc[:,0] 
q2a = - df2a.iloc[:,1] + 1.

df3 = pd.read_csv("VRii.txt", delim_whitespace=True, skiprows=1)

c3 = df3.iloc[:,0] 
q3 = - df3.iloc[:,1] + 0.15

df3a = pd.read_csv("VRiia.txt", delim_whitespace=True, skiprows=1)

c3a = df3a.iloc[:,0] 
q3a = - df3a.iloc[:,1] + 0.15

df4 = pd.read_csv("RIii.txt", delim_whitespace=True, skiprows=1)

c4 = df4.iloc[:,0] 
q4 = - df4.iloc[:,1] + 0.15

df4a = pd.read_csv("RIiia.txt", delim_whitespace=True, skiprows=1)

c4a = df4a.iloc[:,0] 
q4a = - df4a.iloc[:,1] + 0.15

fontsize = 17

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (1, 1, 1)

plt.xlim (0., 1.0)

plt.plot (c1a, q1a, color = 'red',   linewidth = 2, linestyle = 'solid')
plt.plot (c2a, q2a, color = 'green', linewidth = 2, linestyle = 'solid')
plt.plot (c4a, q4a, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.plot (c3a, q3a, color = 'black', linewidth = 2, linestyle = 'solid')

plt.plot ([0.,c1a[0]], [0., q1a[0]], color = 'red',   linewidth = 2, linestyle = 'solid')
plt.plot ([0.,c2a[0]], [0., q2a[0]], color = 'green', linewidth = 2, linestyle = 'solid')
plt.plot ([0.,c3a[0]], [0., q3a[0]], color = 'black', linewidth = 2, linestyle = 'solid')

plt.plot (c1, q1, color = 'red',   linewidth = 2, linestyle = 'dashed', label = "HRi")
plt.plot (c2, q2, color = 'green', linewidth = 2, linestyle = 'dashed', label = "HRii")
plt.plot (c4, q4, color = 'blue',  linewidth = 2, linestyle = 'dashed', label = "RIii")
plt.plot (c3, q3, color = 'black', linewidth = 2, linestyle = 'dashed', label = "VRii")

plt.plot ([0.,c1[0]], [0., q1[0]], color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.plot ([0.,c2[0]], [0., q2[0]], color = 'green', linewidth = 2, linestyle = 'dashed')
plt.plot ([0.,c4[0]], [0., q4[0]], color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot ([0.,c3[0]], [0., q3[0]], color = 'black', linewidth = 2, linestyle = 'dashed')

plt.axhline (0., color = 'black',   linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$c_\beta$',              fontsize = fontsize)
plt.ylabel (r'$(Q_{E} + Q_e)_{crit}$', fontsize = fontsize)
plt.legend (fontsize = fontsize)

plt.tight_layout ()

#plt.show () 
plt.savefig ("Figure6.pdf")
