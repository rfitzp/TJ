# #################################################################
# Script to plot magnetic flux-surface in simple two-filament model
# #################################################################

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

def F (Y):

    f1 = Y*Y
    f2 = (Y+1.)*(Y+1.)

    return 0.5 * math.log (f1) + 0.5 * beta * math.log (f2) - Psi_x

Y_c = root_scalar(F, bracket=[0.01, 1.]).root
X_c = 0.

print ("\nX-point:     (X, Y) = (%10.3e, %10.3e) Psi_x = %10.3e" % (X_x, Y_x, Psi_x))
print ("Start-point: (X, Y) = (%10.3e, %10.3e)\n" % (X_c, Y_c))

yy = input("y: ")
y  = float (yy)

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

print ("q = %11.4e\n" % q)

theta_points = phi_points /q 

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.axis    ('equal')  
plt.plot    (X_points, Y_points, color = "black", linewidth = 2)
plt.axhline (0.,                 color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,                 color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.plot (0., 0.,  'ro')
plt.plot (0., -1., 'ro')
plt.plot (0., Y_x, 'bx')

plt.subplot (2, 2, 2)

plt.plot    (theta_points, phi_points, color = "black", linewidth = 2)
plt.axvline (1.,                       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlim (0., 2.)

plt.ylabel (r'$\phi/\pi$',   fontsize = "15")
plt.xlabel (r'$\theta/\pi$', fontsize = "15")

plt.subplot (2, 2, 3)

plt.plot    (omega_points, phi_points, color = "black", linewidth = 2)
plt.axvline (1.,                       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlim (0., 2.)

plt.ylabel (r'$\phi/\pi$',   fontsize = "15")
plt.xlabel (r'$\omega/\pi$', fontsize = "15")

plt.subplot (2, 2, 4)

plt.plot    (omega_points, theta_points, color = "black", linewidth = 2)
plt.axhline (1.,                         color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (1.,                         color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.ylim (0., 2.)
plt.xlim (0., 2.)

plt.xlabel (r'$\omega/\pi$', fontsize = "15")
plt.ylabel (r'$\theta/\pi$', fontsize = "15")

plt.tight_layout ()
plt.show ()

