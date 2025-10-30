import math
from scipy.optimize import root_scalar
import numpy as np
import scipy.constants as sc
import scipy.special as sp

q0   = 1.01
q95  = 3.5
ntor = 2.

Te = 200.*sc.e
Ti = 200.*sc.e
ne = 3.5e19

Dperp = 1.0

B0 = 5.3
R0 = 6.2
a  = 2.0

Lambda = 23. - math.log ((ne/1.e6)**0.5 /(Te/sc.e))

tau_ee = 6. * 2.**0.5 * math.pi**1.5 * sc.epsilon_0**2.0 * sc.m_e**0.5 * Te**1.5 /Lambda /sc.e**4 /ne
sigmap = 1.96 * ne * sc.e**2 * tau_ee /sc.m_e

tauR = sc.mu_0 * a*a * sigmap
tauP = a*a /Dperp
tauA = (R0 /B0) * (sc.mu_0 * 2.*sc.m_p * ne)**0.5

dbeta = ((5./3.) * 2.*sc.m_p * (Te + Ti))**0.5 /sc.e /B0 /a

Dlayer = tauA**0.5 /tauR**0.25 /tauP**0.25 /dbeta**0.5

def f(x):
    return sp.exp1(q95/x) /sp.exp1(q0/x) - 0.05
    
alpha = root_scalar(f, bracket=[0.1, 100]).root

rhs = ntor * alpha**2 * Dlayer**2

def g(x):
    return x /(-math.log(2.*x)) - rhs

xc = root_scalar(g, bracket=[1.e-4, 0.4999]).root

psic = 1. - sp.exp1 (q0/alpha - math.log (2.*xc)) /sp.exp1 (q0/alpha)

qc = q0 - alpha * math.log (2.*xc)

mc = int (qc * ntor)

print ("\nq0   = %11.4e q95  = %11.4e alpha = %11.4e Lambda = %11.4e"
       % (q0, q95, alpha, Lambda))
print (  "tauA = %11.4e tauP = %11.4e tauR  = %11.4e dbeta  = %11.4e Dlayer = %11.4e"
         % (tauA, tauP, tauR, dbeta, Dlayer))
print (  "rhs  = %11.4e xc   = %11.4e psic  = %11.4e qc     = %11.4e mc     = %3d\n"
         % (rhs, xc, psic, qc, mc))
