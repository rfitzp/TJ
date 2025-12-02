import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("Scannu.txt", delim_whitespace = True, header = None)

m   = df.iloc[:,5]
psi = df.iloc[:,9]
Dr  = df.iloc[:,10]
Di  = df.iloc[:,11]

p3  = []
r3  = []
i3  = []
p4  = []
r4  = []
i4  = []
p5  = []
r5  = []
i5  = []
p6  = []
r6  = []
i6  = []
p7  = []
r7  = []
i7  = []
p8  = []
r8  = []
i8  = []
p9  = []
r9  = []
i9  = []
p10 = []
r10 = []
i10 = []
p11 = []
r11 = []
i11 = []
p12 = []
r12 = []
i12 = []

for mm, p, r, i in zip (m, psi, Dr, Di):
    if mm == 3:
        p3.append(p)
        r3.append(r)
        i3.append(i)
    elif mm == 4:
        p4.append(p)
        r4.append(r)
        i4.append(i)
    elif mm == 5:
        p5.append(p)
        r5.append(r)
        i5.append(i)
    elif mm == 6:
        p6.append(p)
        r6.append(r)
        i6.append(i)
    elif mm == 7:
        p7.append(p)
        r7.append(r)
        i7.append(i)
    elif mm == 8:
        p8.append(p)
        r8.append(r)
        i8.append(i)
    elif mm == 9:
        p9.append(p)
        r9.append(r)
        i9.append(i)
    elif mm == 10:
        p10.append(p)
        r10.append(r)
        i10.append(i)
    elif mm == 11:
        p11.append(p)
        r11.append(r)
        i11.append(i)
    elif mm == 12:
        p12.append(p)
        r12.append(r)
        i12.append(i)   

fontsize = 17
        
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (1, 1, 1)

plt.xlim (0.835, 1.00)
#plt.ylim (0., 250.)

if len(p4) > 0:
    plt.plot (p4,  r4,  color = 'black',   linewidth = 2, linestyle = 'solid', label = '$m=4$')
if len(p5) > 0:
    plt.plot (p5,  r5,  color = 'red',     linewidth = 2, linestyle = 'solid', label = '$m=5$')
if len(p6) > 0:
    plt.plot (p6,  r6,  color = 'green',   linewidth = 2, linestyle = 'solid', label = '$m=6$')
if len(p7) > 0: 
    plt.plot (p7,  r7,  color = 'blue',    linewidth = 2, linestyle = 'solid', label = '$m=7$')
if len(p8) > 0:
    plt.plot (p8,  r8,  color = 'yellow',  linewidth = 2, linestyle = 'solid', label = '$m=8$')
if len(p9) > 0:
    plt.plot (p9,  r9,  color = 'cyan',    linewidth = 2, linestyle = 'solid', label = '$m=9$')    
if len(p10) > 0:
    plt.plot (p10, r10, color = 'magenta', linewidth = 2, linestyle = 'solid', label = '$m=10$')
if len(p11) > 0:
    plt.plot (p11, r11, color = 'brown',   linewidth = 2, linestyle = 'solid', label = '$m=11$')    

if len(p4) > 0:
    plt.plot (p4,  i4,  color = 'black',   linewidth = 2, linestyle = 'dashed')
if len(p5) > 0:
    plt.plot (p5,  i5,  color = 'red',     linewidth = 2, linestyle = 'dashed')
if len(p6) > 0:
    plt.plot (p6,  i6,  color = 'green',   linewidth = 2, linestyle = 'dashed')
if len(p7) > 0: 
    plt.plot (p7,  i7,  color = 'blue',    linewidth = 2, linestyle = 'dashed')
if len(p8) > 0:
    plt.plot (p8,  i8,  color = 'yellow',  linewidth = 2, linestyle = 'dashed')
if len(p9) > 0:
    plt.plot (p9,  i9,  color = 'cyan',    linewidth = 2, linestyle = 'dashed')    
if len(p10) > 0:
    plt.plot (p10, i10, color = 'magenta', linewidth = 2, linestyle = 'dashed')
if len(p11) > 0:
    plt.plot (p11, i11, color = 'brown',   linewidth = 2, linestyle = 'dashed')    
        
plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.ylabel (r'$\Delta_s$', fontsize = fontsize)
plt.xlabel (r"$\Psi_N$",   fontsize = fontsize)

plt.legend (fontsize = fontsize)

plt.tight_layout ()

plt.show ()    

