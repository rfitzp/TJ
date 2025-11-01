# Script to plot safety-factor in simple two-filament model

import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

beta = 0.2
qast = 12.0

def Psi (X, Y):

    f1 = X*X + Y*Y
    f2 = X*X + (Y + 1.) * (Y + 1.)

    return 0.5 * math.log (f1) + 0.5 * beta * math.log (f2)

def PsiX (X, Y):

    f1 = X*X + Y*Y
    f2 = X*X + (Y + 1.) * (Y + 1.)

    return X /f1 + beta * X /f2

def PsiY (X, Y):

    f1 = X*X + Y*Y
    f2 = X*X + (Y + 1.) * (Y + 1.)

    return Y /f1 + beta * (Y + 1.) /f2

def Rhs (t, y):

    X, Y, phi, omega = y

    psiX = PsiX (X, Y)
    psiY = PsiY (X, Y)

    gradPsi = math.sqrt (psiX * psiX + psiY * psiY)
    L2      = X*X + Y*Y

    dXdt     = - psiY /gradPsi
    dYdt     =   psiX /gradPsi
    dphidt   =   qast /gradPsi
    domegadt = (X * psiX + Y * psiY) /L2 /gradPsi

    return [dXdt, dYdt, dphidt, domegadt]

def Stop_Condition (t, y):
    
    X, Y, phi, omega = y

    return omega - 2.*math.pi

Stop_Condition.terminal  = True   
Stop_Condition.direction = +1  

X_x   = 0.
Y_x   = - 1. /(1. + beta)
Psi_x = math.log (beta ** beta /(1. + beta) ** (1. + beta))

def F (X):

    f1 = X*X
    f2 = X*X + 1.

    return 0.5 * math.log (f1) + 0.5 * beta * math.log (f2) - Psi_x

X_c = root_scalar(F, bracket=[0.01, 1.]).root
Y_c = 0.

def Get_q (x):

    X     = x*X_c 
    Y     = 0.0
    phi   = 0.0
    omega = 0.0

    y0 = [X, Y, phi, omega]
    
    sol = solve_ivp(
        Rhs,
        [0., 100.],
        y0,
        method       = 'RK45',
        events       = Stop_Condition,
        dense_output = False,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    qval = sol.y[2][-1] /2./math.pi

    return qval

xx = np.logspace (-1., -1.e-4, 5000, base = 10.0)
qq = []
pp = []
pl = []

n = 0
for x in xx:

    qval = Get_q (x)
    psi  = Psi_x /Psi (x*X_c, Y_c)

    qq.append (qval)
    pp.append (psi)
    pl.append (-math.log10(1.-psi))

    if n%1000 == 0:
        print ("n = %3d" % n)

    n = n + 1

xx1 = np.logspace (1.e-4, 0.1, 5000, base = 10.0)
qq1 = []
pp1 = []
pl1 = []

n = 0
for x in xx1:

    qval = Get_q (x)
    psi  = Psi_x /Psi (x*X_c, Y_c)

    qq1.append (qval)
    pp1.append (psi)
    pl1.append (-math.log10(psi - 1.))

    if n%1000 == 0:
        print ("n = %3d" % n)

    n = n + 1
    
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.plot    (pp,  qq,  color = "red",   linewidth = 2)
plt.plot    (pp1, qq1, color = "blue",  linewidth = 2)
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlim (0.2, 1.5)

plt.xlabel (r'$\Psi$', fontsize = "15")
plt.ylabel (r'$q$',    fontsize = "15")

plt.subplot (2, 1, 2)

plt.plot    (pl,  qq,  color = "red",   linewidth = 2)
plt.plot    (pl1, qq1, color = "blue",  linewidth = 2)
plt.axhline (0.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$-\log_{10}(|1-\Psi|)$', fontsize = "15")
plt.ylabel (r'$q$',                  fontsize = "15")

plt.tight_layout ()
#plt.show ()
plt.savefig ("Figure2.pdf")
