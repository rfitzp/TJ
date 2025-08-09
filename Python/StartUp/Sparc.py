import scipy.constants as sc
import math

Z  = 2.
a  = 0.57
R0 = 1.85
B0 = 12.2
Te = 20.*sc.e*1.e3
ne = 3.1e20

eta = (Z * 15. /1.96/6./2.**0.5/math.pi**1.5) * sc.m_e**0.5 * sc.e**2. * sc.c**4. * sc.mu_0**2. /Te**1.5

tauR = a*a * sc.mu_0 /eta

print ("eta = %11.4e tauR = %11.4e" % (eta, tauR))

Z  = 1.2
a  = 0.21
R0 = 0.67
B0 = 8.0
Te = 3.*sc.e*1.e3
ne = 1.0e20

eta = (Z * 15. /1.96/6./2.**0.5/math.pi**1.5) * sc.m_e**0.5 * sc.e**2. * sc.c**4. * sc.mu_0**2. /Te**1.5

tauR = a*a * sc.mu_0 /eta

print ("eta = %11.4e tauR = %11.4e" % (eta, tauR))
