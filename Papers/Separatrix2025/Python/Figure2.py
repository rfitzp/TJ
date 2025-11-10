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

def Rhs1 (t, y):

    X, Y, phi, omega = y

    psiX = PsiX (X, Y)
    psiY = PsiY (X, Y)

    gradPsi = math.sqrt (psiX * psiX + psiY * psiY)
    L2      = X*X + Y*Y

    dXdt     = 1.
    dYdt     = 0.
    dphidt   = 0.
    domegadt = - Y /L2

    return [dXdt, dYdt, dphidt, domegadt]

def Stop_Condition (t, y):
    
    X, Y, phi, omega = y

    return omega - 2.*math.pi

Stop_Condition.terminal  = True   
Stop_Condition.direction = +1  

X_x   = 0.
Y_x   = - 1. /(1. + beta)
Psi_x = math.log (beta ** beta /(1. + beta) ** (1. + beta))

Y_d = (-1. + Y_x) /2.

def Stop_Condition1 (t, y):
    
    X, Y, phi, omega = y

    return Y - Y_d

Stop_Condition1.terminal  = True   
Stop_Condition1.direction = -1  

def F (Y):

    f1 = Y*Y
    f2 = (Y+1.)*(Y+1.)

    return 0.5 * math.log (f1) + 0.5 * beta * math.log (f2) - Psi_x

Y_c = root_scalar(F, bracket=[0.01, 1.]).root
X_c = 0.

def Get_q (y):
    
    X     = 0.0
    Y     = y*Y_c 
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

def Get_q1 (y):

    X     = 0.0
    Y     = y*Y_c 
    phi   = 0.0
    omega = 0.0

    y0 = [X, Y, phi, omega]

    sol1 = solve_ivp(
        Rhs,
        [0., 100.],
        y0,
        method       = 'RK45',
        events       = Stop_Condition1,
        dense_output = False,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    X_d = sol1.y[0][-1]

    def Stop_Condition2 (t, y):
    
        X, Y, phi, omega = y
        
        return X + X_d
    
    Stop_Condition2.terminal  = True   
    Stop_Condition2.direction = +1
    
    sol2 = solve_ivp(
        Rhs1,
        [0., 100.],
        [sol1.y[0][-1], sol1.y[1][-1], sol1.y[2][-1], sol1.y[3][-1]],
        method       = 'RK45',
        events       = Stop_Condition2,
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    sol3 = solve_ivp(
        Rhs,
        [0., 100.],
        [sol2.y[0][-1], sol2.y[1][-1], sol2.y[2][-1], sol2.y[3][-1]],
        method       = 'RK45',
        events       = Stop_Condition,
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )

    qval = sol3.y[2][-1] /2./math.pi

    return qval

yy = np.logspace (-1., -1.e-5, 5000, base = 10.0)
qq = []
pp = []
pl = []

n = 0
for y in yy:

    qval = Get_q (y)
    psi  = Psi_x /Psi (X_c, y*Y_c)

    qq.append (qval)
    pp.append (psi)
    pl.append (-math.log10(1.-psi))

    if n%1000 == 0:
        print ("n = %3d" % n)

    n = n + 1

yy1 = np.logspace (1.e-4, 0.1, 5000, base = 10.0)
qq1 = []
pp1 = []
pl1 = []

n = 0
for y in yy1:

    qval = Get_q (y)
    psi  = Psi_x /Psi (X_c, y*Y_c)

    qq1.append (qval)
    pp1.append (psi)
    pl1.append (-math.log10(psi - 1.))

    if n%1000 == 0:
        print ("n = %3d" % n)

    n = n + 1

qq2 = []
pp2 = []
pl2 = []

n = 0
for y in yy1:

    qval = Get_q1 (y)
    psi  = Psi_x /Psi (X_c, y*Y_c)

    qq2.append (qval)
    pp2.append (psi)
    pl2.append (-math.log10(psi - 1.))

    if n%1000 == 0:
        print ("n = %3d" % n)

    n = n + 1    
    
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (2, 1, 1)

plt.plot    (pp,  qq,  color = "red",   linewidth = 2)
plt.plot    (pp1, qq1, color = "blue",  linewidth = 2)
plt.plot    (pp2, qq2, color = "green", linewidth = 2)
plt.axvline (1.,       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlim (0.2, 1.5)
plt.ylim (0.,  7.5)

plt.xlabel (r'$\Psi$', fontsize = "20")
plt.ylabel (r'$q$',    fontsize = "20")

plt.subplot (2, 1, 2)

plt.plot    (pl,  qq,  color = "red",   linewidth = 2)
plt.plot    (pl1, qq1, color = "blue",  linewidth = 2)
plt.plot    (pl2, qq2, color = "green", linewidth = 2)

plt.xlim (0., 3.25)
plt.ylim (0., 7.5)

plt.xlabel (r'$-\log_{10}(|\Psi-1|)$', fontsize = "20")
plt.ylabel (r'$q$',                    fontsize = "20")

plt.tight_layout ()
#plt.show ()
plt.savefig ("Figure2.pdf")
