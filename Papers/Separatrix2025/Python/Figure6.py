# Script to plot outer flux coordinates in simple two-filament model

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

def GetPointsOuter (y):

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
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    tmax1     = sol1.t_events[0][0]
    t_points1 = np.linspace (0., tmax1, 1000)
    
    X_points1     = sol1.sol (t_points1)[0]
    Y_points1     = sol1.sol (t_points1)[1]
    phi_points1   = sol1.sol (t_points1)[2] 
    omega_points1 = sol1.sol (t_points1)[3] 

    X_d = X_points1[-1]

    def Stop_Condition2 (t, y):
    
        X, Y, phi, omega = y
        
        return X + X_d
    
    Stop_Condition2.terminal  = True   
    Stop_Condition2.direction = +1
    
    sol2 = solve_ivp(
        Rhs1,
        [0., 100.],
        [X_points1[-1], Y_points1[-1], phi_points1[-1], omega_points1[-1]],
        method       = 'RK45',
        events       = Stop_Condition2,
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )
    
    tmax2     = sol2.t_events[0][0]
    t_points2 = np.linspace (0., tmax2, 100)
    
    X_points2     = sol2.sol (t_points2)[0]
    Y_points2     = sol2.sol (t_points2)[1]
    phi_points2   = sol2.sol (t_points2)[2] 
    omega_points2 = sol2.sol (t_points2)[3]
    
    sol3 = solve_ivp(
        Rhs,
        [0., 100.],
        [X_points2[-1], Y_points2[-1], phi_points2[-1], omega_points2[-1]],
        method       = 'RK45',
        events       = Stop_Condition,
        dense_output = True,
        rtol         = 1e-8,
        atol         = 1e-10
    )

    tmax3     = sol3.t_events[0][0]
    t_points3 = np.linspace (0., tmax3, 1000)
    
    X_points3     = sol3.sol (t_points3)[0]
    Y_points3     = sol3.sol (t_points3)[1]
    phi_points3   = sol3.sol (t_points3)[2]
    omega_points3 = sol3.sol (t_points3)[3]
    
    q = phi_points3[-1] /2./math.pi
    
    theta_points1 = phi_points1 /q
    theta_points2 = phi_points2 /q
    theta_points3 = phi_points3 /q

    t_points2 += t_points1[-1]
    t_points3 += t_points2[-1]

    X_points     = np.concatenate ((np.array(X_points1),     np.array(X_points2),     np.array(X_points3)))
    Y_points     = np.concatenate ((np.array(Y_points1),     np.array(Y_points2),     np.array(Y_points3)))
    phi_points   = np.concatenate ((np.array(phi_points1),   np.array(phi_points2),   np.array(phi_points3)))
    omega_points = np.concatenate ((np.array(omega_points1), np.array(omega_points2), np.array(omega_points3)))
    theta_points = np.concatenate ((np.array(theta_points1), np.array(theta_points2), np.array(theta_points3)))
    t_points     = np.concatenate ((np.array(t_points1),     np.array(t_points2),     np.array(t_points3)))

    return q, X_points, Y_points, phi_points, omega_points, theta_points, t_points

npsi   = 1000
ntheta = 80
nt     = 1000

yy = np.logspace (1.e-5, 0.1,           npsi, base = 10.0)
tt = np.linspace (0., 2.*math.pi-1.e-8, ntheta)

Xvals = np.empty((npsi, nt))
Yvals = np.empty((npsi, nt))
Xval1 = np.empty((npsi, ntheta))
Yval1 = np.empty((npsi, ntheta))

n = 0
for y in yy:
    
    q, X_points, Y_points, phi_points, omega_points, theta_points, t_points = GetPointsOuter (y)

    tx = np.linspace (0., t_points[-1], nt)
    
    Xint = interp1d (t_points, X_points, kind = 'linear')
    Yint = interp1d (t_points, Y_points, kind = 'linear')
    
    Xvals[n,:] = Xint (tx)
    Yvals[n,:] = Yint (tx)
    
    Xint1 = interp1d (theta_points, X_points, kind = 'linear')
    Yint1 = interp1d (theta_points, Y_points, kind = 'linear')

    Xnew = Xint1 (tt)
    Ynew = Yint1 (tt)

    Xval1[n,:] = Xnew
    Yval1[n,:] = Ynew

    n = n + 1

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 1, 1)

plt.axis ('equal')

for n in range (0, npsi, 50):
    plt.plot (Xvals[n,:], Yvals[n,:], color = "red", linewidth = 1.)
plt.plot (Xvals[-1,:], Yvals[-1,:],   color = "red", linewidth = 1.)
for n in range (ntheta):
    plt.plot (Xval1[:,n], Yval1[:,n], color = "blue", linewidth = 1.)
    
plt.plot (0., 0.,  'ko')

plt.plot ((-0.3, 0.3), (Y_d, Y_d),  color = 'black', linewidth = 3.)


plt.xlabel (r'$X$', fontsize = "20")
plt.ylabel (r'$Y$', fontsize = "20")

plt.tight_layout ()
#plt.show ()
plt.savefig ("Figure6.pdf")
    
