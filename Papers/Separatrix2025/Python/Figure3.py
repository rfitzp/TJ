# Script to plot inner flux coordinates in simple two-filament model

import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
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

def F (Y):

    f1 = Y*Y
    f2 = (Y+1.)*(Y+1.)

    return 0.5 * math.log (f1) + 0.5 * beta * math.log (f2) - Psi_x

Y_c = root_scalar(F, bracket=[0.01, 1.]).root
X_c = 0.

npsi   = 1000
ntheta = 80
nt     = 1000

yy = np.logspace (-0.5, -1.e-5,         npsi, base = 10.0)
tt = np.linspace (0., 2.*math.pi-1.e-8, ntheta)

Xvals = np.empty((npsi, nt))
Yvals = np.empty((npsi, nt))
Xval1 = np.empty((npsi, ntheta))
Yval1 = np.empty((npsi, ntheta))

n = 0
for y in yy:
    
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
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    tmax     = sol.t_events[0][0]
    t_points = np.linspace (0., tmax, nt)
    
    X_points   = sol.sol (t_points)[0]
    Y_points   = sol.sol (t_points)[1]
    phi_points = sol.sol (t_points)[2]

    Xvals[n,:] = X_points
    Yvals[n,:] = Y_points

    q = phi_points[-1] /2./math.pi

    theta_points = phi_points /q 

    Xint = interp1d (theta_points, X_points, kind = 'linear')
    Yint = interp1d (theta_points, Y_points, kind = 'linear')

    Xnew = Xint (tt)
    Ynew = Yint (tt)

    Xval1[n,:] = Xnew
    Yval1[n,:] = Ynew

    n = n + 1

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (1, 1, 1)

plt.axis ('equal')

for n in range (0, npsi, 10):
    plt.plot (Xvals[n,:], Yvals[n,:], color = "red", linewidth = 1.)
plt.plot (Xvals[-1,:], Yvals[-1,:],   color = "red", linewidth = 1.)
for n in range (ntheta):
    plt.plot (Xval1[:,n], Yval1[:,n], color = "blue", linewidth = 1.)
    
#plt.axhline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')
#plt.axvline (0., color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (0., 0.,  'ko')
#plt.plot (0., -1., 'ro')

plt.tight_layout ()
#plt.show ()
plt.savefig("Figure3.pdf")
    
