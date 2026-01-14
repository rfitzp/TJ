# Script to plot magnetic flux-surface in simple two-filament model with divertor

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

def GetPointsInner (y):

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
    t_points = np.linspace (0., tmax, 1000)
    
    X_points     = sol.sol (t_points)[0]
    Y_points     = sol.sol (t_points)[1]
    phi_points   = sol.sol (t_points)[2] /math.pi
    omega_points = sol.sol (t_points)[3] /math.pi
    
    q = phi_points[-1] /2.

    theta_points = phi_points /q

    return q, X_points, Y_points, phi_points, omega_points, theta_points

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
    phi_points1   = sol1.sol (t_points1)[2] /math.pi
    omega_points1 = sol1.sol (t_points1)[3] /math.pi

    X_d = X_points1[-1]

    def Stop_Condition2 (t, y):
    
        X, Y, phi, omega = y
        
        return X + X_d
    
    Stop_Condition2.terminal  = True   
    Stop_Condition2.direction = +1
    
    sol2 = solve_ivp(
        Rhs1,
        [0., 100.],
        [X_points1[-1], Y_points1[-1], math.pi*phi_points1[-1], math.pi*omega_points1[-1]],
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
    phi_points2   = sol2.sol (t_points2)[2] /math.pi
    omega_points2 = sol2.sol (t_points2)[3] /math.pi
    
    sol3 = solve_ivp(
        Rhs,
        [0., 100.],
        [X_points2[-1], Y_points2[-1], math.pi*phi_points2[-1], math.pi*omega_points2[-1]],
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
    phi_points3   = sol3.sol (t_points3)[2] /math.pi
    omega_points3 = sol3.sol (t_points3)[3] /math.pi
    
    q = phi_points3[-1] /2.
    
    theta_points1 = phi_points1 /q
    theta_points2 = phi_points2 /q
    theta_points3 = phi_points3 /q

    X_points     = np.concatenate ((np.array(X_points1),     np.array(X_points2),     np.array(X_points3)))
    Y_points     = np.concatenate ((np.array(Y_points1),     np.array(Y_points2),     np.array(Y_points3)))
    phi_points   = np.concatenate ((np.array(phi_points1),   np.array(phi_points2),   np.array(phi_points3)))
    omega_points = np.concatenate ((np.array(omega_points1), np.array(omega_points2), np.array(omega_points3)))
    theta_points = np.concatenate ((np.array(theta_points1), np.array(theta_points2), np.array(theta_points3)))

    return q, X_points, Y_points, phi_points, omega_points, theta_points

q1, X_points1, Y_points1, phi_points1, omega_points1, theta_points1 = GetPointsInner (0.94533)
q2, X_points2, Y_points2, phi_points2, omega_points2, theta_points2 = GetPointsOuter (1.0)
q3, X_points3, Y_points3, phi_points3, omega_points3, theta_points3 = GetPointsOuter (1.047)

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 1, 1)

plt.plot ((-0.3, 0.3), (Y_d, Y_d),  color = 'black', linewidth = 3.)

plt.axis    ('equal')  
plt.plot    (X_points1, Y_points1, color = "red",   linewidth = 2)
plt.plot    (X_points3, Y_points3, color = "green", linewidth = 2)
plt.plot    (X_points2, Y_points2, color = "blue",  linewidth = 2)
plt.axhline (0.,                   color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,                   color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (0., 0.,  'ko')

plt.xlabel (r'$X$', fontsize = "20")
plt.ylabel (r'$Y$', fontsize = "20")

plt.tight_layout ()
#plt.show ()
plt.savefig ("Figure5.pdf")
