# Ideal.py

# Plots ideal response parameter versus PSI for nu scan

import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("../../Outputs/TearX/Scannu.txt", delim_whitespace = True, header = None)

psi = df.iloc[:,6]
m   = df.iloc[:,5]
id  = np.log10(np.asarray(df.iloc[:,12]))

p3  = []
i3  = []
p4  = []
i4  = []
p5  = []
i5  = []
p6  = []
i6  = []
p7  = []
i7  = []
p8  = []
i8  = []
p9  = []
i9  = []
p10 = []
i10 = []
p11 = []
i11 = []

for p, mm, i in zip (psi, m, id):
    if mm == 3:
        p3.append(p)
        i3.append(i)
    elif mm == 4:
        p4.append(p)
        i4.append(i)
    elif mm == 5:
        p5.append(p)
        i5.append(i)
    elif mm == 6:
        p6.append(p)
        i6.append(i)
    elif mm == 7:
        p7.append(p)
        i7.append(i)
    elif mm == 8:
        p8.append(p)
        i8.append(i)
    elif mm == 9:
        p9.append(p)
        i9.append(i)
    elif mm == 10:
        p10.append(p)
        i10.append(i)
    elif mm == 11:
        p11.append(p)
        i11.append(i)   

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'TEARX Code: nu Scan')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0.8, 1.00)

if len(p4) > 0:
    plt.plot (p4, i4,   color = 'black',   linewidth = 2, linestyle = 'solid', label = '$m=4$')
if len(p5) > 0:
    plt.plot (p5, i5,   color = 'red',     linewidth = 2, linestyle = 'solid', label = '$m=5$')
if len(p6) > 0:
    plt.plot (p6, i6,   color = 'green',   linewidth = 2, linestyle = 'solid', label = '$m=6$')
if len(p7) > 0: 
    plt.plot (p7, i7,   color = 'blue',    linewidth = 2, linestyle = 'solid', label = '$m=7$')
if len(p8) > 0:
    plt.plot (p8, i8,   color = 'yellow',  linewidth = 2, linestyle = 'solid', label = '$m=8$')
if len(p9) > 0:
    plt.plot (p9, i9,   color = 'cyan',    linewidth = 2, linestyle = 'solid', label = '$m=9$')    
if len(p10) > 0:
    plt.plot (p10, i10, color = 'magenta', linewidth = 2, linestyle = 'solid', label = '$m=10$')
if len(p11) > 0:
    plt.plot (p11, i11, color = 'brown',   linewidth = 2, linestyle = 'solid', label = '$m=11$')    
    
plt.axhline (0.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',                         fontsize = "15")
plt.ylabel (r"$\log_{10}[|\Delta_s|/(-E_{ss})]$", fontsize = "15")
plt.legend (fontsize = '15')

plt.tight_layout ()

plt.show ()    
