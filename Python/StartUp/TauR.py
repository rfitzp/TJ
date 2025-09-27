# ################################################
# Script to simulate current ramp in ohmic tokamak
# ################################################

import scipy.constants as sc
import numpy as np
import math
import matplotlib.pyplot as plt

Z  = 2.
Ln = 15.
Te = 7.e3*sc.e

eta = (Z * Ln /1.96 /6. /2.**0.5 /math.pi**1.5) * (sc.m_e**0.5 * sc.e**2 * sc.c**4 * sc.mu_0**2) /Te**1.5

# ###
# JET
# ###
a = 0.96

tauR1 = a*a * sc.mu_0 /eta

# #####
# SPARC
# #####
a = 0.57

tauR2 = a*a * sc.mu_0 /eta

# ####
# ITER
# ####
a  = 2.0

tauR3 = a*a * sc.mu_0 /eta

# ####
# DEMO
# ####
a = 2.93

tauR4 = a*a * sc.mu_0 /eta

print ("\ntau_R: JET - %11.4e SPARC - %11.4e ITER - %11.4e DEMO - %11.4e\n" % (tauR1, tauR2, tauR3, tauR4))
