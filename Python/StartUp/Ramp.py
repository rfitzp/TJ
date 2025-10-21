# ################################################
# Script to simulate current ramp in ohmic tokamak
# ################################################

import scipy.constants as sc
import numpy as np
import math
import matplotlib.pyplot as plt

T_0 = ( (15. /1.96 /6./2.**0.5 /math.pi**1.5) * (sc.m_e**0.5 * sc.c**4 * sc.e**2 /1.e20) )**0.4 /sc.e /1.e3

beta_p = T_0 * sc.e * 1.e3 * 1.e20 * sc.mu_0

t_0 = 1./beta_p

I_0 = 2.*math.pi /sc.mu_0 /1.e6

P_0 = 4.*math.pi**2 * 1.e20 * T_0 * sc.e * 1.e3 /1.e6

E_0 = beta_p

E_c = sc.e**3 * 1.e20 * 15. /4./math.pi /sc.epsilon_0**2 /sc.m_e /sc.c**2

E_D = E_c * sc.m_e * sc.c**2 /T_0/sc.e/1.e3

print ("\nT_0 = %11.4e keV beta_p = %11.4e t_0 = %11.4e s I_0 = %11.4e kappa MA  P_0 = %11.4e MW"
       % (T_0, beta_p, t_0, I_0, P_0))
print ("\nE_c = %11.4e V/m E_D = %11.4e V/m E_c_hat = %11.4e E_D_hat = %11.4e"
       % (E_c, E_D, E_c/E_0, E_D/E_0))

# .........
# alpha = 0
# .........
qa    = 3.304
Tramp = 2.078
Eramp = 2.205
tramp = 0.2137

chi   = 1.0
Z     = 2.0

# ####
# ITER
# ####
B0    = 5.3
R0    = 6.2
a     = 2.0
n20   = 1.0
kappa = 1.85
delta = 0.0
chi   = 1.0

gs  = 1. + kappa*kappa * (1. + 2.*delta**2. - 1.2*delta**3.) /2.
gp  = (gs/kappa)**2.
gnc = 4.3 - 0.6 * R0/a
gnc = 1.
Z0  = gnc * Z

Bt = (a/R0) * B0 /qa

nG = (5./math.pi) * gs * Bt /a

#for i in range(5):
bp = beta_p            * Z0**(+0.4) * n20**(0.6)  * chi**(-0.4) * Bt**(-1.2)
T0 = T_0               * Z0**(+0.4) * n20**(-0.4) * chi**(-0.4) * Bt**(+0.8)
t0 = t_0     * a*a     * Z0**(-0.4) * n20**(-0.6) * chi**(-0.6) * Bt**(+1.2)
I0 = 5.00    * a                                                * Bt        * gs
P0 = P_0     * R0      * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(+0.8)* gp
E0 = beta_p  * a**(-1) * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(-0.2)
Ec = E_c     *                        n20
ED = E_D               * Z0**(-0.4) * n20**(+1.4) * chi**(0.4)  * Bt**(-0.8)

tauE = 0.0562 * I0**0.83 * B0**0.15 * (P0*Eramp)**(-0.69) * n20**0.41 * R0**1.97 * kappa**0.78 * (a/R0)**0.58
chi  = a*a /tauE

print ("\nITER: gs   = %11.4e gp = %11.4e gnc = %11.4e" % (gs, gp, gnc))
print ("\nITER: T0   = %11.4e keV t0 = %11.4e s I0 = %11.4e MA E0 = %11.4e V/m Ec = %11.4e V/m ED = %11.4e V/m"
       % (T0, t0, I0, E0, Ec, ED))
print ("\nITER: Tc   = %11.4e keV t_ramp = %11.4e s E = %11.4e V/m P = %11.4e (MW) beta_p = %11.4e nG = %11.4e 10^20 m^-3"
       % (T0*Tramp, t0*tramp, E0*Eramp, P0*Eramp, bp, nG))
print ("\nITER: tauE = %11.4e s   chi    = %11.4e m^2/s" % (tauE, chi))

# ##########
# Simulation
# ##########
tpre = 0.05

tt = np.linspace (0.01, 1., 5000)

xx  = []
dd  = []
IIp = []
qqa = []
TTc = []
PPc = []
EE  = []
EEc = []
EED = []
nnG = []

for t in tt:

    xx.append (t*t0*tramp)

    if (t < tpre):
        dd .append (tpre**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa*tpre/t)
        TTc.append (T0*Tramp * tpre**0.4    * (t/tpre)**0.8)
        PPc.append (P0*Eramp * tpre**0.4    * (t/tpre)**0.8)
        EE .append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2))
        EEc.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /Ec)
        EED.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /ED /(T0*Tramp * tpre**0.4 * (t/tpre)**0.8))
        nnG.append ((n20/nG)  * (t/tpre)**(-1.))
    else:
        dd .append (t**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa)
        TTc.append (T0*Tramp * t**0.4)
        PPc.append (P0*Eramp * t**0.4)
        EE .append (E0*Eramp * t**(-0.6))
        EEc.append (E0*Eramp * t**(-0.6) /Ec)
        EED.append (E0*Eramp * t**(-0.6) /ED /(T0*Tramp * t**0.4))
        nnG.append (n20/nG) 

font = 15
fig  = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 2, 1)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, dd, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.ylabel (r'$a (m)$', fontsize = font)

plt.subplot (3, 2, 2)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, IIp, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = "$I_p(MA)$")
plt.plot    (xx, PPc, color = 'red',   linewidth = 2, linestyle = 'solid',  label = "$P(MW)$")
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 3)

plt.xlim (0., xx[-1])

plt.plot    (xx, qqa, color = 'blue',  linewidth = 2, linestyle = 'solid', label = "$q_a$")
plt.plot    (xx, nnG, color = 'red',   linewidth = 2, linestyle = 'solid', label = "$n_e/n_G$")
plt.axhline (0.,      color = 'black', linewidth = 0, linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',  fontsize = font)
#plt.ylabel (r'$q_a$',    fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 4)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, TTc, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',     fontsize = font)
plt.ylabel (r'$T_c (keV)$', fontsize = font)

plt.subplot (3, 2, 5)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EE, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',          fontsize = font)
plt.ylabel (r'${\cal E} (V/m)$', fontsize = font)

plt.subplot (3, 2, 6)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EEc, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_c$")
plt.plot    (xx, EED, color = 'red',   linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_D$")
plt.axhline (0.,      color = 'black', linewidth = 1,  linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()

#plt.show () 
plt.savefig("Figure6.jpg", dpi=300)
plt.close()

# #####
# SPARC
# #####

B0    = 12.2
R0    = 1.85
a     = 0.57
n20   = 2.0
kappa = 1.97
delta = 0.0
chi   = 1.0

gs  = 1. + kappa*kappa * (1. + 2.*delta**2. - 1.2*delta**3.) /2.
gp  = (gs/kappa)**2.
#gnc = 4.3 - 0.6 * R0/a
gnc = 1.
Z0  = gnc * Z

Bt = (a/R0) * B0 /qa

nG = (5./math.pi) * gs * Bt /a

#for i in range(5):
bp = beta_p            * Z0**(+0.4) * n20**(0.6)  * chi**(-0.4) * Bt**(-1.2)
T0 = T_0               * Z0**(+0.4) * n20**(-0.4) * chi**(-0.4) * Bt**(+0.8)
t0 = t_0     * a*a     * Z0**(-0.4) * n20**(-0.6) * chi**(-0.6) * Bt**(+1.2)
I0 = 5.00    * a                                                * Bt        * gs
P0 = P_0     * R0      * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(+0.8)* gp
E0 = beta_p  * a**(-1) * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(-0.2)
Ec = E_c     *                        n20
ED = E_D               * Z0**(-0.4) * n20**(+1.4) * chi**(0.4)  * Bt**(-0.8)

tauE = 0.0562 * I0**0.83 * B0**0.15 * (P0*Eramp)**(-0.69) * n20**0.41 * R0**1.97 * kappa**0.78 * (a/R0)**0.58
chi  = a*a /tauE

print ("\nSPARC: gs = %11.4e gp = %11.4e gnc = %11.4e" % (gs, gp, gnc))
print ("\nSPARC: T0 = %11.4e keV t0 = %11.4e s I0 = %11.4e MA E0 = %11.4e V/m Ec = %11.4e V/m ED = %11.4e V/m"
       % (T0, t0, I0, E0, Ec, ED))
print ("\nSPARC: Tc = %11.4e keV t_ramp = %11.4e s E = %11.4e V/m P = %11.4e MW beta_p = %11.4e nG = %11.4e 10^20 m^-3"
       % (T0*Tramp, t0*tramp, E0*Eramp, P0*Eramp, bp, nG))
print ("\nSPARC: tauE = %11.4e s   chi  = %11.4e m^2/s" % (tauE, chi))

# ##########
# Simulation
# ##########
tpre = 0.05

tt = np.linspace (0.01, 1., 5000)

xx  = []
dd  = []
IIp = []
qqa = []
TTc = []
PPc = []
EE  = []
EEc = []
EED = []
nnG = []

for t in tt:

    xx.append (t*t0*tramp)

    if (t < tpre):
        dd .append (tpre**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa*tpre/t)
        TTc.append (T0*Tramp * tpre**0.4    * (t/tpre)**0.8)
        PPc.append (P0*Eramp * tpre**0.4    * (t/tpre)**0.8)
        EE .append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2))
        EEc.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /Ec)
        EED.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /ED /(T0*Tramp * tpre**0.4 * (t/tpre)**0.8))
        nnG.append ((n20/nG)  * (t/tpre)**(-1.))
    else:
        dd .append (t**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa)
        TTc.append (T0*Tramp * t**0.4)
        PPc.append (P0*Eramp * t**0.4)
        EE .append (E0*Eramp * t**(-0.6))
        EEc.append (E0*Eramp * t**(-0.6) /Ec)
        EED.append (E0*Eramp * t**(-0.6) /ED /(T0*Tramp * t**0.4))
        nnG.append (n20/nG) 

font = 15
fig  = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 2, 1)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, dd, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.ylabel (r'$a (m)$', fontsize = font)

plt.subplot (3, 2, 2)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, IIp, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = "$I_p(MA)$")
plt.plot    (xx, PPc, color = 'red',   linewidth = 2, linestyle = 'solid',  label = "$P(MW)$")
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',    fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 3)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, qqa, color = 'blue',  linewidth = 2, linestyle = 'solid', label = "$q_a$")
plt.plot    (xx, nnG, color = 'red',   linewidth = 2, linestyle = 'solid', label = "$n_e/n_G$")
plt.axhline (0.,      color = 'black', linewidth = 0, linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',  fontsize = font)
#plt.ylabel (r'$q_a$',    fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 4)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, TTc, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',     fontsize = font)
plt.ylabel (r'$T_c (keV)$', fontsize = font)

plt.subplot (3, 2, 5)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EE, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',          fontsize = font)
plt.ylabel (r'${\cal E} (V/m)$', fontsize = font)

plt.subplot (3, 2, 6)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EEc, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_c$")
plt.plot    (xx, EED, color = 'red',   linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_D$")
plt.axhline (0.,      color = 'black', linewidth = 1,  linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()

#plt.show ()
plt.savefig("Figure5.jpg", dpi=300)
plt.close()

# ###
# JET
# ###

B0    = 3.45
R0    = 2.96
a     = 0.96
n20   = 0.3
kappa = 1.8
delta = 0.
chi   = 1.0

gs  = 1. + kappa*kappa * (1. + 2.*delta**2. - 1.2*delta**3.) /2.
gp  = (gs/kappa)**2.
#gnc = 4.3 - 0.6 * R0/a
gnc = 1.
Z0  = gnc * Z

Bt = (a/R0) * B0 /qa

nG = (5./math.pi) * gs * Bt /a

#for i in range(5):
bp = beta_p            * Z0**(+0.4) * n20**(0.6)  * chi**(-0.4) * Bt**(-1.2)
T0 = T_0               * Z0**(+0.4) * n20**(-0.4) * chi**(-0.4) * Bt**(+0.8)
t0 = t_0     * a*a     * Z0**(-0.4) * n20**(-0.6) * chi**(-0.6) * Bt**(+1.2)
I0 = 5.00    * a                                                * Bt        * gs
P0 = P_0     * R0      * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(+0.8)* gp
E0 = beta_p  * a**(-1) * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(-0.2)
Ec = E_c     *                        n20
ED = E_D               * Z0**(-0.4) * n20**(+1.4) * chi**(0.4)  * Bt**(-0.8)

tauE = 0.0562 * I0**0.83 * B0**0.15 * (P0*Eramp)**(-0.69) * n20**0.41 * R0**1.97 * kappa**0.78 * (a/R0)**0.58
chi  = a*a /tauE
    
print ("\nJET: gs   = %11.4e gp = %11.4e gnc = %11.4e" % (gs, gp, gnc))
print ("\nJET: T0   = %11.4e keV t0 = %11.4e s I0 = %11.4e MA E0 = %11.4e V/m Ec = %11.4e V/m ED = %11.4e V/m"
       % (T0, t0, I0, E0, Ec, ED))
print ("\nJET: Tc   = %11.4e keV t_ramp = %11.4e s E = %11.4e V/m P = %11.4e MW beta_p = %11.4e nG = %11.4e 10^20 m^-3"
       % (T0*Tramp, t0*tramp, E0*Eramp, P0*Eramp, bp, nG))
print ("\nJET: tauE = %11.4e s  chi   = %11.4e m^2/s" % (tauE, chi))

# ##########
# Simulation
# ##########
tpre = 0.05

tt = np.linspace (0.01, 1., 5000)

xx  = []
dd  = []
IIp = []
qqa = []
TTc = []
PPc = []
EE  = []
EEc = []
EED = []
nnG = []

for t in tt:

    xx.append (t*t0*tramp)

    if (t < tpre):
        dd .append (tpre**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa*tpre/t)
        TTc.append (T0*Tramp * tpre**0.4    * (t/tpre)**0.8)
        PPc.append (P0*Eramp * tpre**0.4    * (t/tpre)**0.8)
        EE .append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2))
        EEc.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /Ec)
        EED.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /ED /(T0*Tramp * tpre**0.4 * (t/tpre)**0.8))
        nnG.append ((n20/nG)  * (t/tpre)**(-1.))
    else:
        dd .append (t**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa)
        TTc.append (T0*Tramp * t**0.4)
        PPc.append (P0*Eramp * t**0.4)
        EE .append (E0*Eramp * t**(-0.6))
        EEc.append (E0*Eramp * t**(-0.6) /Ec)
        EED.append (E0*Eramp * t**(-0.6) /ED /(T0*Tramp * t**0.4))
        nnG.append (n20/nG) 

font = 15
fig  = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 2, 1)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, dd, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.ylabel (r'$a (m)$', fontsize = font)

plt.subplot (3, 2, 2)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, IIp, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = "$I_p(MA)$")
plt.plot    (xx, PPc, color = 'red',   linewidth = 2, linestyle = 'solid',  label = "$P(MW)$")
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',    fontsize = font)
#plt.ylabel (r'$I_p (MA)$', fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 3)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, qqa, color = 'blue',  linewidth = 2, linestyle = 'solid', label = "$q_a$")
plt.plot    (xx, nnG, color = 'red',   linewidth = 2, linestyle = 'solid', label = "$n_e/n_G$")
plt.axhline (0.,      color = 'black', linewidth = 0, linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',  fontsize = font)
#plt.ylabel (r'$q_a$',    fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 4)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, TTc, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',     fontsize = font)
plt.ylabel (r'$T_c (keV)$', fontsize = font)

plt.subplot (3, 2, 5)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EE, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',          fontsize = font)
plt.ylabel (r'${\cal E} (V/m)$', fontsize = font)

plt.subplot (3, 2, 6)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EEc, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_c$")
plt.plot    (xx, EED, color = 'red',   linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_D$")
plt.axhline (0.,      color = 'black', linewidth = 1,  linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()

#plt.show ()
plt.savefig("Figure4.jpg", dpi=300)
plt.close()

# ###
# DEMO
# ###

B0    = 5.86
R0    = 9.07
a     = 2.93
n20   = 0.88
kappa = 1.65
delta = 0.
chi   = 1.0

gs  = 1. + kappa*kappa * (1. + 2.*delta**2. - 1.2*delta**3.) /2.
gp  = (gs/kappa)**2.
#gnc = 4.3 - 0.6 * R0/a
gnc = 1.
Z0  = gnc * Z

Bt = (a/R0) * B0 /qa

nG = (5./math.pi) * gs * Bt /a

#for i in range(5):
bp = beta_p            * Z0**(+0.4) * n20**(0.6)  * chi**(-0.4) * Bt**(-1.2)
T0 = T_0               * Z0**(+0.4) * n20**(-0.4) * chi**(-0.4) * Bt**(+0.8)
t0 = t_0     * a*a     * Z0**(-0.4) * n20**(-0.6) * chi**(-0.6) * Bt**(+1.2)
I0 = 5.00    * a                                                * Bt        * gs
P0 = P_0     * R0      * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(+0.8)* gp
E0 = beta_p  * a**(-1) * Z0**(+0.4) * n20**(+0.6) * chi**(+0.6) * Bt**(-0.2)
Ec = E_c     *                        n20
ED = E_D               * Z0**(-0.4) * n20**(+1.4) * chi**(0.4)  * Bt**(-0.8)

tauE = 0.0562 * I0**0.83 * B0**0.15 * (P0*Eramp)**(-0.69) * n20**0.41 * R0**1.97 * kappa**0.78 * (a/R0)**0.58
chi  = a*a /tauE

print ("\nDEMO: gs   = %11.4e gp = %11.4e gnc = %11.4e" % (gs, gp, gnc))
print ("\nDEMO: T0   = %11.4e keV t0 = %11.4e s I0 = %11.4e MA E0 = %11.4e V/m Ec = %11.4e V/m ED = %11.4e V/m"
       % (T0, t0, I0, E0, Ec, ED))
print ("\nDEMO: Tc   = %11.4e keV t_ramp = %11.4e s E = %11.4e V/m P = %11.4e MW beta_p = %11.4e nG = %11.4e 10^20 m^-3"
       % (T0*Tramp, t0*tramp, E0*Eramp, P0*Eramp, bp, nG))
print ("\nDEMO: tauE = %11.4e s   chi    = %11.4e m^2/s" % (tauE, chi))

# ##########
# Simulation
# ##########
tpre = 0.05

tt = np.linspace (0.01, 1., 5000)

xx  = []
dd  = []
IIp = []
qqa = []
TTc = []
PPc = []
EE  = []
EEc = []
EED = []
nnG = []

for t in tt:

    xx.append (t*t0*tramp)

    if (t < tpre):
        dd .append (tpre**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa*tpre/t)
        TTc.append (T0*Tramp * tpre**0.4    * (t/tpre)**0.8)
        PPc.append (P0*Eramp * tpre**0.4    * (t/tpre)**0.8)
        EE .append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2))
        EEc.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /Ec)
        EED.append (E0*Eramp * tpre**(-0.6) * (t/tpre)**(-0.2) /ED /(T0*Tramp * tpre**0.4 * (t/tpre)**0.8))
        nnG.append ((n20/nG)  * (t/tpre)**(-1.))
    else:
        dd .append (t**0.5 * a)
        IIp.append (t * I0)
        qqa.append (qa)
        TTc.append (T0*Tramp * t**0.4)
        PPc.append (P0*Eramp * t**0.4)
        EE .append (E0*Eramp * t**(-0.6))
        EEc.append (E0*Eramp * t**(-0.6) /Ec)
        EED.append (E0*Eramp * t**(-0.6) /ED /(T0*Tramp * t**0.4))
        nnG.append (n20/nG) 

font = 15
fig  = plt.figure (figsize = (8.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (3, 2, 1)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, dd, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.ylabel (r'$a (m)$', fontsize = font)

plt.subplot (3, 2, 2)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, IIp, color = 'blue',  linewidth = 2, linestyle = 'solid',  label = "$I_p(MA)$")
plt.plot    (xx, PPc, color = 'red',   linewidth = 2, linestyle = 'solid',  label = "$P(MW)$")
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',    fontsize = font)
#plt.ylabel (r'$I_p (MA)$', fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 3)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, qqa, color = 'blue',  linewidth = 2, linestyle = 'solid', label = "$q_a$")
plt.plot    (xx, nnG, color = 'red',   linewidth = 2, linestyle = 'solid', label = "$n_e/n_G$")
plt.axhline (0.,      color = 'black', linewidth = 0, linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',  fontsize = font)
#plt.ylabel (r'$q_a$',    fontsize = font)
plt.legend (fontsize = font)

plt.subplot (3, 2, 4)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, TTc, color = 'blue',  linewidth = 2, linestyle = 'solid')
plt.axhline (0.,      color = 'black', linewidth = 1, linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',     fontsize = font)
plt.ylabel (r'$T_c (keV)$', fontsize = font)

plt.subplot (3, 2, 5)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EE, color = 'blue',  linewidth = 2,  linestyle = 'solid')
plt.axhline (0.,     color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$',          fontsize = font)
plt.ylabel (r'${\cal E} (V/m)$', fontsize = font)

plt.subplot (3, 2, 6)

plt.xlim (0., xx[-1])
 
plt.plot    (xx, EEc, color = 'blue',  linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_c$")
plt.plot    (xx, EED, color = 'red',   linewidth = 2,  linestyle = 'solid', label = "${\cal E}/E_D$")
plt.axhline (0.,      color = 'black', linewidth = 1,  linestyle = 'dotted')
plt.axhline (1.,      color = 'black', linewidth = 1,  linestyle = 'dotted')

plt.axvline (tpre*t0*tramp, color = 'black', linewidth = 1, linestyle = 'dotted')

plt.xlabel (r'$t (s)$', fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()

#plt.show ()
plt.savefig("Figure7.jpg", dpi=300)
plt.close()
