import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

df = pd.read_csv ("Scannu.txt", delim_whitespace = True, header = None)

id  = df.iloc[:,9]
m   = df.iloc[:,5]
psi = df.iloc[:,4]

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
    elif mm == 12:
        p12.append(p)
        i12.append(i)   

fontsize = 17
        
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize) 

plt.subplot (2, 1, 1)

plt.ylim (0.835, 1.00)
plt.xlim (2., 9.)

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
    
#plt.axhline (1.,    color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$q_{95}$', fontsize = fontsize)
plt.ylabel (r"$\Psi_N$", fontsize = fontsize)

plt.legend (fontsize = 16)

psi = df.iloc[:,9]
idd = df.iloc[:,8]

xp3  = []
xi3  = []
xp4  = []
xi4  = []
xp5  = []
xi5  = []
xp6  = []
xi6  = []
xp7  = []
xi7  = []
xp8  = []
xi8  = []
xp9  = []
xi9  = []
xp10 = []
xi10 = []
xp11 = []
xi11 = []
xp12 = []
xi12 = []

for p, mm, i in zip (psi, m, idd):
    if mm == 3:
        xp3.append(p)
        xi3.append(i)
    elif mm == 4:
        xp4.append(p)
        xi4.append(i)
    elif mm == 5:
        xp5.append(p)
        xi5.append(i)
    elif mm == 6:
        xp6.append(p)
        xi6.append(i)
    elif mm == 7:
        xp7.append(p)
        xi7.append(i)
    elif mm == 8:
        xp8.append(p)
        xi8.append(i)
    elif mm == 9:
        xp9.append(p)
        xi9.append(i)
    elif mm == 10:
        xp10.append(p)
        xi10.append(i)   
    elif mm == 11:
        xp11.append(p)
        xi11.append(i)
    elif mm == 12:
        xp12.append(p)
        xi12.append(i)

plt.subplot (2, 1, 2)

plt.xlim (0.835, 1.00)

if len(p4) > 0:
    plt.plot (xp4, xi4,   color = 'black',   linewidth = 2, linestyle = 'solid', label = '$m=4$')
if len(p5) > 0:
    plt.plot (xp5, xi5,   color = 'red',     linewidth = 2, linestyle = 'solid', label = '$m=5$')
if len(p6) > 0:
    plt.plot (xp6, xi6,   color = 'green',   linewidth = 2, linestyle = 'solid', label = '$m=6$')
if len(p7) > 0: 
    plt.plot (xp7, xi7,   color = 'blue',    linewidth = 2, linestyle = 'solid', label = '$m=7$')
if len(p8) > 0:
    plt.plot (xp8, xi8,   color = 'yellow',  linewidth = 2, linestyle = 'solid', label = '$m=8$')
if len(p9) > 0:
    plt.plot (xp9, xi9,   color = 'cyan',    linewidth = 2, linestyle = 'solid', label = '$m=9$')    
if len(p10) > 0:
    plt.plot (xp10, xi10, color = 'magenta', linewidth = 2, linestyle = 'solid', label = '$m=10$')
if len(p11) > 0:
    plt.plot (xp11, xi11, color = 'brown',   linewidth = 2, linestyle = 'solid', label = '$m=11$')
    
plt.ticklabel_format (style = 'sci', axis = 'y', scilimits = (0, 0))
plt.xlabel (r'$\Psi_N$', fontsize = fontsize)
plt.ylabel (r"$E_{ss}$", fontsize = fontsize)
#plt.legend (fontsize = fontsize)

plt.tight_layout ()

plt.show ()    
#plt.savefig ("Figure9.pdf")
