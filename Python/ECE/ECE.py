# ECE.py

# #############################################
# Script to calculate ECE convolution functions
# #############################################

import math
import scipy.special as sp
import scipy.constants as sc
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

sqpi = math.sqrt(math.pi)

# ##################
# Physics parameters
# ##################
Te = 5.     # Electron temperature at resonance (keV)
ne = 1.e20  # Electron number density at resonance (m^-3)
B0 = 5.3    # On-axis magnetic field-strength (T)
R0 = 6.2    # Major radius of magnetic axis (m)
Rw = 6.5    # Major radius of resonance (m)

# ##################
# Derived parameters
# ##################
wc0  = sc.e * B0 /sc.m_e                                 # On-axis cyclotron frequency
tau0 = wc0 * R0 /sc.c                                    # Optical depth parameter
vtc2 = Te*1.e3 * sc.e /sc.c/sc.c /sc.m_e                 # (v_t /c)^2 at resonance
wp   = math.sqrt (ne * sc.e*sc.e /sc.epsilon_0 /sc.m_e)  # Electron plasma frequency at resonance
w1   =     (wc0 /wp) * (R0 /Rw)                          # Ratio of 1st harmonic resonant frequency to plasma frequency
w2   = 2.* (wc0 /wp) * (R0 /Rw)                          # Ratio of 2nd harmonic resonant frequency to plasma frequency
tw   = Te*1.e3 * sc.e /sc.m_e /sc.c/sc.c                 # Relativistic temperature parameter

print ('w_c0 = %11.4e tau0 = %11.4e (vt/c)^2 = %11.4e w1_hat^2 = %11.4e w2_hat^2 = %11.4e tw = %11.4e' % (wc0, tau0, vtc2, w1*w1, w2*w2, tw))

# ###################
# Absorption function
# ###################
def Get_F72 (z):

    if (z > 0):
        sqz = math.sqrt(z)
    else:
        sqz = math.sqrt(-z) * 1j
        
    f72 = (8./15.) * (0.75 - 0.5*z + z*z - sqpi * z*z*sqz * sp.wofz(sqz*1j))

    return f72.real, f72.imag

# #################################
# First harmonic resonance function
# #################################
def Get_z1(wc):

    return (1. - wc/w1) /vtc2

# ##################################
# Second harmonic resonance function
# ##################################
def Get_z2(wc):

    return (1. - 2.*wc/w2) /vtc2

# #######################################
# Refractive index of 1st harmonic O-mode
# #######################################
def Get_NperpO(wc):

    z1         = Get_z1(wc)
    F72r, F72i = Get_F72(z1)

    wpc2  = 1. /wc/wc
    denom = 1. + 0.5 * wpc2 * F72r
    N2    = (1. - wpc2) /denom

    return math.sqrt (N2)

# #######################################
# Refractive index of 2nd harmonic X-mode
# #######################################
def Get_NperpX(wc):

    z2         = Get_z2(wc)
    F72r, F72i = Get_F72(z2)

    wpc2  = 1. /wc/wc

    NC2 = 1. - (wpc2 /3.) * (1. - 0.25 * wpc2) /(1. - wpc2 /3.)
    a   = - 0.5 * wpc2 * F72r / (1. - wpc2 /3.)
    b   = - 2. * (1. - wpc2 /6.) * a
    N2  = NC2 * (1. - (b + a * NC2))
    
    return math.sqrt (N2)

# #############################################
# Absorption coefficient of 1st harmonic O-mode
# #############################################
def Get_alpha1O(wc):

    z1         = Get_z1(wc)
    F72r, F72i = Get_F72(z1)

    wpc2   = 1. /wc/wc
    denom  = 1. + 0.5 * wpc2 * F72r
    NperpO = Get_NperpO(wc)
 
    return 0.5 * NperpO * wpc2 * (-F72i) /denom

# #############################################
# Absorption coefficient of 2nd harmonic X-mode
# #############################################
def Get_alpha2X(wc):

    z2         = Get_z2(wc)
    F72r, F72i = Get_F72(z2)

    wpc2   = 1. /wc/wc
    NperpX = Get_NperpX(wc)
    N2     = NperpX * NperpX
    a2     = 0.5 * wpc2 * (1. + 3. * N2 * F72r) /(3. - wpc2 * (1. + 1.5 * N2 * F72r))
    A2     = (1. + a2) **2
    denom  = 1. + 0.5 * wpc2 * A2 * F72r

    return NperpX * wpc2 * A2 * (-F72i) /denom

# ###############################################
# Integrand for 1st harmonic O-mode optical depth
# ###############################################
def Get_f1O(wc):

    return tau0 * Get_alpha1O(wc) /wc

# ###################################
# First harmonic O-mode optical depth
# ###################################
def Get_tau1O(wc):

    if w1 < wc:
        res, err = quad (Get_f1O, w1, wc)
        return res
    else:
        return 0.

# ###############################################    
# Integrand for 2nd harmonic X-mode optical depth
# ###############################################
def Get_f2X(wc):

    return tau0 * Get_alpha2X(wc) /wc

# ####################################
# Second harmonic X-mode optical depth
# ####################################
def Get_tau2X(wc):

    if w2/2. < wc:       
        res, err = quad (Get_f2X, w2/2., wc)
        return res
    else:
        return 0.

# #############################    
# Calculate absorption function
# #############################
zmax = 20.
znum = 1000

zz = np.linspace(-zmax, zmax, znum)

f72r = []
f72i = []
for z in zz:
    F72r, F72i = Get_F72(z)
    f72r.append (F72r)
    f72i.append (F72i)

# ##################################################    
# Calculate 1st harmonic O-mode convolution function
# ##################################################
wcOmin = 0.98*w1
wcOmax = 1.1*w1
wcOnum = 1000

wwcO = np.linspace(wcOmin, wcOmax, wcOnum)

alphaO = []
tauO   = []
GO     = []
RO     = []
HO     = []
for wcO in wwcO:

    alpha = Get_alpha1O(wcO)
    tau   = Get_tau1O(wcO)
    G     = tau0 * alpha * math.exp(-tau) /wcO
    R     = R0 * wc0 /wp /wcO
    H     = (wc0 /wp) * G * R0 /R/R
    
    alphaO.append (alpha)
    tauO.append (tau)
    GO.append (G)
    RO.append (R)
    HO.append (H)

# ##################################################    
# Calculate 2nd harmonic X-mode convolution function
# ##################################################
wcXmin = 0.98*w2/2.
wcXmax = 1.1*w2/2.
wcXnum = 1000

wwcX = np.linspace(wcXmin, wcXmax, wcXnum)

alphaX = []
tauX   = []
GX     = []
RX     = []
HX     = []
for wcX in wwcX:

    alpha = Get_alpha2X(wcX)
    tau   = Get_tau2X(wcX)
    G     = tau0 * alpha * math.exp(-tau) /wcX
    R     = R0 * wc0 /wp /wcX
    H     = (wc0 /wp) * G * R0 /R/R

    alphaX.append (alpha)
    tauX.append (tau)
    GX.append (G)
    RX.append (R)
    HX.append (H)

# #################################################################    
# Find center and widths of O-mode and X-mode convolution functions
# #################################################################
HOmax = max(HO)
HXmax = max(HX)

index = HO.index(max(HO))
ROmax = RO[index]
index = HX.index(max(HX))
RXmax = RX[index]

sigmaO = 1./math.sqrt(2.*math.pi) /HOmax
sigmaX = 1./math.sqrt(2.*math.pi) /HXmax

# ########################################################
# Gaussian fit to 1st harmonic O-mode convolution function
# ########################################################
HOm = []
for R in RO:

    h0 = math.exp (- (R - ROmax) * (R - ROmax) /2./sigmaO/sigmaO) /math.sqrt(2.*math.pi) /sigmaO
    
    HOm.append (h0)

# ########################################################    
# Gaussian fit to 2nd harmonic X-mode convolution function
# ########################################################
HXm = []
for R in RX:
    
    h0 = math.exp (- (R - RXmax) * (R - RXmax) /2./sigmaX/sigmaX) /math.sqrt(2.*math.pi) /sigmaX

    HXm.append (h0)    
 
print ("sigmaO = %11.4e  (Rw-ROmax)/sigmaX = %11.4e sigmaX = %11.4e (Rw-RXmax)/sigmaO = %11.4e" %
       (sigmaO, (Rw-ROmax)/sigmaO, sigmaX, (Rw-RXmax)/sigmaX))

# #####
# Plots
# #####
"""
font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.xlim (-zmax, zmax)
 
plt.plot    (zz, f72r, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = r'$Re(F_{7/2})$')
plt.plot    (zz, f72i, color = 'red',   linewidth = 2,  linestyle = 'solid', label = r'$Im(F_{7/2})$')
plt.axhline (0.,       color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (0.,       color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$z$',fontsize = font)
plt.legend (fontsize = font)

plt.savefig("F72.pdf")
"""

"""
font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font)

plt.subplot (2, 2, 1)

plt.xlim (wcOmin, wcOmax)
 
plt.plot    (wwcO, alphaO, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,           color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w1,           color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\hat{\alpha}_1^{\,(O)}$', fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$',      fontsize = font)

plt.subplot (2, 2, 2)

plt.xlim (wcOmin, wcOmax)
 
plt.plot    (wwcO, tauO, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,         color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w1,         color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\tau_1^{\,(O)}$',    fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$', fontsize = font)

plt.subplot (2, 2, 3)

plt.xlim (wcOmin, wcOmax)
 
plt.plot    (wwcO, GO, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,       color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w1,       color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\omega_p\,H_1^{\,(0)}$', fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$',     fontsize = font)

plt.subplot (2, 2, 4)

plt.xlim (R0 * wc0 /wp /wcOmax, R0 * wc0 /wp /wcOmin)
 
plt.plot    (RO, HO,           color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.plot    (RO, HOm,          color = 'red',   linewidth = 2,  linestyle = 'dotted')
plt.axhline (0.,               color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (R0 * wc0 /wp /w1, color = 'black', linewidth = 1., linestyle = 'dashed')
plt.axvline (R0,               color = 'black', linewidth = 1., linestyle = 'dotted')

plt.ylabel (r'$F_1^{\,(0)}$', fontsize = font)
plt.xlabel (r'$R$',           fontsize = font)

plt.tight_layout ()

plt.savefig("O.pdf")
"""

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font)

plt.subplot (2, 2, 1)

plt.xlim (wcXmin, wcXmax)
 
plt.plot    (wwcX, alphaX, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,           color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w2/2.,        color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\hat{\alpha}_2^{\,(X)}$', fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$',      fontsize = font)

plt.subplot (2, 2, 2)

plt.xlim (wcXmin, wcXmax)
 
plt.plot    (wwcX, tauX, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,         color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w2/2.,      color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\tau_2^{\,(X)}$',    fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$', fontsize = font)

plt.subplot (2, 2, 3)

plt.xlim (wcXmin, wcXmax)
 
plt.plot    (wwcX, GX, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,       color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (w2/2.,    color = 'black', linewidth = 1., linestyle = 'dashed')

plt.ylabel (r'$\omega_p\,H_2^{\,(X)}$', fontsize = font)
plt.xlabel (r'$\omega_c/\omega_p$',     fontsize = font)

plt.subplot (2, 2, 4)

plt.xlim (R0 * wc0 /wp /wcXmax, R0 * wc0 /wp /wcXmin)
 
plt.plot    (RX, HX,           color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.plot    (RX, HXm,          color = 'red',   linewidth = 2,  linestyle = 'dotted')
plt.axhline (0.,               color = 'black', linewidth = 1., linestyle = 'dotted')
plt.axvline (R0 * wc0 /wp /w1, color = 'black', linewidth = 1., linestyle = 'dashed')
plt.axvline (R0,               color = 'black', linewidth = 1., linestyle = 'dotted')

plt.ylabel (r'$F_2^{\,(X)}$', fontsize = font)
plt.xlabel (r'$R$',           fontsize = font)

plt.tight_layout ()

plt.savefig("X.pdf")    
#plt.show()

