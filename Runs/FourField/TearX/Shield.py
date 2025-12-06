import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("Scannu.txt", delim_whitespace = True, header = None)

m  = df.iloc[:,5]
r  = df.iloc[:,6]
E  = df.iloc[:,8]
Dr = df.iloc[:,10]
Di = df.iloc[:,11]

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
p12 = []
i12 = []

for mm, rr, e, dr, di in zip (m, r, E, Dr, Di):

    i = ((1. + dr/(-e)) * (1. + dr/(-e)) + di*di/e/e)**0.5 - 1.
    
    if mm == 3:
        p3.append(rr)
        i3.append(i)
    elif mm == 4:
        p4.append(rr)
        i4.append(i)
    elif mm == 5:
        p5.append(rr)
        i5.append(i)
    elif mm == 6:
        p6.append(rr)
        i6.append(i)
    elif mm == 7:
        p7.append(rr)
        i7.append(i)
    elif mm == 8:
        p8.append(rr)
        i8.append(i)
    elif mm == 9:
        p9.append(rr)
        i9.append(i)
    elif mm == 10:
        p10.append(rr)
        i10.append(i)   
    elif mm == 11:
        p11.append(rr)
        i11.append(i)
    elif mm == 12:
        p12.append(rr)
        i12.append(i)

fontsize = 17        

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0.8, 1.00)
plt.ylim (0., 240.)

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
    
plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$\hat{r}$', fontsize = fontsize)
plt.ylabel (r"$\Sigma$",  fontsize = fontsize)

#plt.legend (fontsize = fontsize)

plt.tight_layout ()

#plt.show ()    
plt.savefig("Figure11.pdf")
