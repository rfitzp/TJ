import math
import scipy.special as sp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import numpy as np

q0  = 1.
q95 = 3.5

def f(x):
    return sp.exp1(q95/x) /sp.exp1(q0/x) - 0.05
    
alpha = root_scalar(f, bracket=[0.1, 100]).root

r95 = (1. - math.exp(-(q95 - q0)/alpha))**0.5

def q(x):
    return q0 + (q95 - q0) * math.log(1.-x*x) /math.log(1.-r95*r95)

def psi(x):
    return 1. - sp.exp1(q0/alpha - math.log(1.-x**2)) /sp.exp1(q0/alpha)

def s(x):
    return 2.*alpha*x*x/(1.-x*x)/q(x)

start = 0.8
stop  = 0.999

rr = np.linspace(start,stop,1000)

qq = []
pp = []
ss = []                                   

for r in rr:
    qq.append(q(r))
    pp.append(psi(r))
    ss.append(s(r))                                 

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 1, 1)

plt.xlim (start, stop)

plt.plot (rr, qq, color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (q95, color = 'red', linewidth = 1.0, linestyle = 'dotted')
plt.axvline (r95, color = 'red', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$q$',       fontsize = "15")

plt.subplot (3, 1, 2)

plt.xlim (start, stop)

plt.plot (rr, pp, color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axhline (0.95, color = 'red', linewidth = 1.0, linestyle = 'dotted')
plt.axvline (r95,  color = 'red', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$',    fontsize = "15")
plt.ylabel (r'$\hat{\psi}$', fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (start, stop)

plt.plot (rr, ss, color = 'blue', linewidth = 2, linestyle = 'solid')

plt.axvline (r95, color = 'red', linewidth = 1.0, linestyle = 'dotted')

plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$s$',       fontsize = "15")                                   

plt.tight_layout ()

plt.show ()    
