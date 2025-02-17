# Solve.py

# Checks analytic expressions for aLL, aLS, aSL, and aSS
# against solutions of diffenential equation and also
# against Bishop paper expressions

import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma
import matplotlib.pyplot as plt

def aLL(nu):

    g1 = gamma(1.  + nu/2.)
    g2 = gamma(1.5 + nu)
    g3 = gamma(0.5 + nu/2.)
    g4 = gamma(1.  + nu)

    return 2.**(1.+nu) * g1 * g2 /(1. + nu) /g3 /g4

def aSL(nu):

    g1 = gamma(1.  + nu/2.)
    g2 = gamma(1.5 + nu)
    g3 = gamma(0.5 + nu/2.)
    g4 = gamma(1.  + nu)

    return  - 2.**nu * nu * g3 * g2 /g1 /g4

def aLS(nu):

    g1 = gamma(0.5 - nu/2.)
    g2 = gamma(0.5 - nu)
    g3 = gamma(1.  - nu)
    g4 = gamma(1.  - nu/2.)

    return  - nu * g1 * g2 /2.**(1.+nu) /g3 /g4

def aSS(nu):

    g1 = gamma(0.5 - nu/2.)
    g2 = gamma(0.5 - nu)
    g3 = gamma(1.  - nu)
    g4 = gamma(1.  - nu/2.)

    return (1. + nu) * g4 * g2 /2.**nu /g1 /g3

def Differential (x, Y):
    
    u, v = Y
    dudx = v
    dvdx = nu*(1.+nu) * u /(1. + x*x)
    
    return [dudx, dvdx]

nu = float (input ("nu ? "))
xmax = 100.

u =       xmax**(-nu)    -             (nu*(1.+nu) /2./(3.+2.*nu)) * xmax**(-nu-2.)
v = -nu * xmax**(-nu-1.) -  (-nu-2.) * (nu*(1.+nu) /2./(3.+2.*nu)) * xmax**(-nu-3.)

Y0 = [u, v]

x_span = (xmax, 0.)
x_eval = np.linspace(xmax, 0., 100)

sol = solve_ivp (Differential, x_span, Y0, t_eval=x_eval)

print ('a_LL = (%11.4e, %11.4e)' % (sol.y[0][-1], aLL(nu)))
print ('a_SL = (%11.4e, %11.4e)' % (sol.y[1][-1], aSL(nu)))

u1 =           xmax**(nu+1.) -            (nu*(1.+nu) /2./(1.-2.*nu)) * xmax**(nu-1.)
v1 = (nu+1.) * xmax**(nu   ) -  (nu-1.) * (nu*(1.+nu) /2./(1.-2.*nu)) * xmax**(nu-2.)

Y01 = [u1, v1]

x_span1 = (xmax, 0.)
x_eval1 = np.linspace(xmax, 0., 100)

sol1 = solve_ivp (Differential, x_span1, Y01, t_eval=x_eval1)

print ('a_LS = (%11.4e, %11.4e)' % (sol1.y[0][-1], aLS(nu)))
print ('a_SS = (%11.4e, %11.4e)' % (sol1.y[1][-1], aSS(nu)))

F1 = 2.**(2.*nu+1.) * gamma(1.-nu) * gamma(0.5+nu) * gamma(1.5+nu) * math.cos(math.pi*nu) /math.pi /gamma(2.+nu)
F2 = (1.+nu) * gamma(1.-nu/2.) * gamma(0.5+nu/2.) /gamma(1.+nu/2.) /gamma(0.5-nu/2.)

tan = math.tan(math.pi*nu/2.)

exp1 = - tan * F2
exp2 =   F2/F1

fitzp1  = aSL(nu) /aLL(nu)
fitzp2  = aSS(nu) /aLL(nu)

print ()
print ('aSL/aLL = (%11.4e, %11.4e)' % (exp1, fitzp1))
print ('aSS/aLL = (%11.4e, %11.4e)' % (exp2, fitzp2))

"""
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 1, 1)

plt.xlim (0., xmax)
 
plt.plot    (sol.t, sol.y[0], color = 'red',   linewidth = 2,   linestyle = 'solid', label = r'$\psi$')
plt.plot    (sol.t, sol.y[1], color = 'green', linewidth = 2,   linestyle = 'solid', label = r'$d\psi/dx$')
plt.axhline (0.,    color = 'black', linewidth = 1.,  linestyle = 'dotted')

plt.xlabel (r'$X$',    fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., xmax)
 
plt.plot    (sol1.t, sol1.y[0], color = 'red',   linewidth = 2,   linestyle = 'solid', label = r'$\psi$')
plt.plot    (sol1.t, sol1.y[1], color = 'green', linewidth = 2,   linestyle = 'solid', label = r'$d\psi/dx$')
plt.axhline (0.,    color = 'black', linewidth = 1.,  linestyle = 'dotted')

plt.xlabel (r'$X$',    fontsize = "15")
plt.legend (fontsize = "15")


plt.tight_layout ()

plt.show ()    
"""
