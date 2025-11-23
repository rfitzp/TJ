import math
import scipy.special as sp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import numpy as np

q0   = 1.01
q95  = 3.5
q105 = 4.0
ntor = 1.

def fm(x):
    return sp.exp1(q95 /x) /sp.exp1(q0 /x) - 0.05
    
alpham = root_scalar(fm, bracket = [0.1, 100]).root

r95 = (1. - math.exp(- (q95 - q0) /alpham))**0.5

def fp(x):
    return sp.exp1(q105 /x) /x - math.exp (q0 /alpham) * sp.exp1(q95 /alpham) /alpham

alphap = root_scalar(fp, bracket = [0.1, 10.]).root

r105 = (1. + math.exp (- q105 /alphap))**0.5

print ("alpha_m = %10.3e alpha_p = %10.3e r_95 = %10.3e r_105 = %10.3e" % (alpham, alphap, r95, r105))

def qm(x):
    return q0 - alpham * math.log(1. - x*x)

def qp(x):
    return - alphap * math.log(x*x - 1.)

def psim(x):
    return 1. - sp.exp1(q0 /alpham - math.log(1. - x**2)) /sp.exp1(q0 /alpham)

def psip(x):
    return 1. + (alpham /alphap) * math.exp (- q0 /alpham) * sp.exp1(- math.log(x**2 - 1.)) /sp.exp1(q0 /alpham)

rm = np.linspace(0.9, 0.999999, 10000)

qqm = []
ppm = []

for r in rm:
    qqm.append(qm(r))
    ppm.append(psim(r))

rp = np.linspace(1.0001, 1.1, 10000)

qqp = []
ppp = []

for r in rp:
    qqp.append(qp(r))
    ppp.append(psip(r))
    
mstart = int (q0     /ntor) + 1    
mstop  = int (qqm[-1]/ntor)

rmm = []
pmm = []

for mpol in range (mstart, mstop, 1):

    qs = float (mpol) /ntor
    
    rmm.append ((1. - math.exp((q0 - qs) /alpham))**0.5)
    pmm.append (1. - sp.exp1(qs /alpham) /sp.exp1(q0 /alpham))

mstart = int (qqp[0] /ntor) 
mstop  = int (qqp[-1]/ntor) 

rmp = []
pmp = []

for mpol in range (mstop, mstart, 1):

    qs = float (mpol) /ntor
    rx = (1. + math.exp (- qs /alphap))**0.5
    px = psip (rx)
    rmp.append (rx)
    pmp.append (px)

fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.xlim (0.95, 1.05)

for r in rmm:
    plt.axvline (r, color = 'black', linewidth = 0.5, linestyle = 'solid')
for r in rmp:
    plt.axvline (r, color = 'black', linewidth = 0.5, linestyle = 'solid')

plt.plot (rm, qqm, color = 'red', linewidth = 2, linestyle = 'solid')
plt.plot (rp, qqp, color = 'red', linewidth = 2, linestyle = 'solid')

plt.ylim (3., 18.)

#plt.xticks([0.95,0.975,1.0,1.025,1.05])
plt.xlabel (r'$\hat{r}$', fontsize = "15")
plt.ylabel (r'$q$',       fontsize = "15")

plt.tight_layout ()

plt.show ()    
plt.savefig ("Figure7.pdf")
