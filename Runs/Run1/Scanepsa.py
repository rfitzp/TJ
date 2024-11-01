import math
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import glob
import os
from scipy.optimize import curve_fit

def model_function(x, a, b, c):
    return a * x**2 + b * x + c

def calc_err(xx, yy, a, b, c):
    sum = 0.
    for x, y in zip(xx, yy):
        res = abs(1. - (a*x*x + b*x + c)/y)
        sum = sum + res
    return sum /float(len(xx))   

directory_path = "."
eq_files = glob.glob(f"{directory_path}/Equilibrium*.nc")
eq_files.sort()
tj_files = glob.glob(f"{directory_path}/TJ*.nc")
tj_files.sort()

eps = []

for eq_file in eq_files:
    ds   = nc.Dataset(eq_file)
    para = ds['InputParameters']
    epsa = float(para[5])
    eps.append(epsa)

E11r = []
E11i = []
E12r = []
E12i = []
E21r = []
E21i = []
E22r = []
E22i = []
    
for tj_file in tj_files:
    ds    = nc.Dataset(tj_file)
    ematr = ds['Emat_r']
    emati = ds['Emat_i']
    E11r.append(float(ematr[0,0]))
    E11i.append(float(emati[0,0]))
    E12r.append(float(ematr[0,1]))
    E12i.append(float(emati[0,1]))
    E21r.append(float(ematr[1,0]))
    E21i.append(-float(emati[1,0]))
    E22r.append(float(ematr[1,1]))
    E22i.append(float(emati[1,1]))

file = open ("fit.txt", "w")    
params, covariance = curve_fit (model_function, eps, E11r)
a, b, c = params
print (f"E11_r = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E11r, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E11i)
a, b, c = params
print (f"E11_i = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E11i, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E12r)
a, b, c = params
print (f"E12_r = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E12r, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E21r)
a, b, c = params
print (f"E21_r = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E21r, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E12i)
a, b, c = params
print (f"E12_i = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E21i, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E21i)
a, b, c = params
print (f"E21_i = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E21i, a, b, c)), file=file)    

params, covariance = curve_fit (model_function, eps, E22r)
a, b, c = params
print (f"E22_r = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E22r, a, b, c)), file=file)

params, covariance = curve_fit (model_function, eps, E22i)
a, b, c = params
print (f"E22_i = %+10.3e %+10.3e*epsa %+10.3e*epsa*epsa: err = %8.1e" % (c, b, a, calc_err (eps, E22i, a, b, c)), file=file)

file.close()

fig = plt.figure (figsize = (12.0, 8.0))
fig.canvas.manager.set_window_title (r'Aspect-Ratio Scan')
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 2, 1)

plt.xlim (0., 0.3)

plt.plot (eps, E11r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)

plt.xlabel (r'$\epsilon_a$', fontsize = "15")
plt.ylabel (r"$Re(E_{11})$", fontsize = "15")

plt.subplot (3, 2, 2)

plt.xlim (0., 0.3)

plt.plot (eps, E11i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon_a$', fontsize = "15")
plt.ylabel (r"$Im(E_{11})$", fontsize = "15")

plt.subplot (3, 2, 3)

plt.xlim (0., 0.3)

plt.plot (eps, E12r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
plt.plot (eps, E21r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)

plt.xlabel (r'$\epsilon_a$',             fontsize = "15")
plt.ylabel (r"$Re(E_{12}), Re(E_{21})$", fontsize = "15")

plt.subplot (3, 2, 4)

plt.xlim (0., 0.3)

plt.plot (eps, E12i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)
plt.plot (eps, E21i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 'x', fillstyle = 'none', markersize = 5)

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon_a$',              fontsize = "15")
plt.ylabel (r"$Im(E_{12}), -Im(E_{21})$", fontsize = "15")

plt.subplot (3, 2, 5)

plt.xlim (0., 0.3)

plt.plot (eps, E22r, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)

plt.xlabel (r'$\epsilon_a$', fontsize = "15")
plt.ylabel (r"$Re(E_{22})$", fontsize = "15")

plt.subplot (3, 2, 6)

plt.xlim (0., 0.3)

plt.plot (eps, E22i, color = 'green', linewidth = 1, linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5)

plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon_a$', fontsize = "15")
plt.ylabel (r"$Im(E_{22})$", fontsize = "15")

plt.tight_layout ()

plt.show ()    
#plt.savefig("Test3.pdf")

