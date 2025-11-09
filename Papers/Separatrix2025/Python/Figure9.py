# ##########################################################################
# Script to plot profile data close to magnetic separatrix in improved model
# ##########################################################################

import math
import numpy as np
import scipy.constants as sc
import scipy.special as sp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# #############################
# JET 84600 NF 61 016001 (2021)
# #############################
q0   = 1.01
q95  = 3.5
q105 = 4.0

B0 = 3.45
R0 = 2.96
a  = 1.25

Wd    = 0.02
Teped = 400.
Tesep = 100.
Tesol = 10.

Tiped = 400.
Tisep = 100.
Tisol = 10.

neped = 5e19
nesep = 2e19
nesol = 0.5e19

Z = 10.
M = 2.

chiperp = 1.
chiphi  = 1.

omegaEfac = 0.

# .........
# Mode data
# .........
ntor = 1.

# .........
# Plot data
# .........
rmin = 0.994
rmax = 1.008

# .................................
# Defininition of mtanh() functions
# .................................
def GetfD (fped, fsep, fsol, fW):

    return fW * math.atanh ((2.*fsep - fped - fsol) /(fped - fsol))

def mtanh (psi, fped, fsep, fsol, fW):

    fD = GetfD (fped, fsep, fsol, fW)

    return fped + ((fsol - fped) /2.) * (1. + math.tanh ((psi - 1. - fD) /fW))

def mtanhp (psi, fped, fsep, fsol, fW):

    fD = GetfD (fped, fsep, fsol, fW)

    return ((fsol - fped) /2.) * (1./math.cosh ((psi - 1. - fD) / fW))**2 /fW

# ..............................
# Calculate alpha_minus and r_95
# ..............................
def fm(x):

    return sp.exp1 (q95 /x) /sp.exp1 (q0 /x) - 0.05
    
alpham = root_scalar (fm, bracket = [0.1, 100]).root

r95 = (1. - math.exp(- (q95 - q0) /alpham))**0.5

# ..............................
# Calculate alpha_plus and r_105
# ..............................
def fp(x):

    return sp.exp1 (q105 /x) /x - math.exp (q0 /alpham) * sp.exp1 (q95 /alpham) /alpham

alphap = root_scalar (fp, bracket = [0.1, 10.]).root

r105 = (1. + math.exp (- q105 /alphap))**0.5

# .....................
# Safety-factor profile
# .....................
def qm(x):

    return q0 - alpham * math.log (1. - x*x)

def qp(x):

    return - alphap * math.log (x*x - 1.)

def q(x):

    if x < 1.:
        return qm(x)
    else:
        return qp(x)

# ......................
# Magnetic shear profile
# ......................
def qmp(x):

    return   alpham * 2. * x*x /(1. - x*x)

def qpp(x):

    return - alphap * 2. * x*x /(x*x - 1.)

def shear(x):

    qval = q(x)
    
    if x < 1:
        return qmp(x) /qval
    else:
        return qpp(x) /qval
    
# .......................
# Normalized flux profile
# .......................
def Psim(x):

    return 1. - sp.exp1 (q0 /alpham - math.log(1. - x**2)) /sp.exp1 (q0 /alpham)

def Psip(x):

    return 1. + (alpham /alphap) * math.exp (- q0 /alpham) * sp.exp1 (- math.log(x**2 - 1.)) /sp.exp1 (q0 /alpham)

def Psi(x):

    if x < 1.:
        return Psim(x)
    else:
        return Psip(x)

# ............................
# Electron temperature profile
# ............................
def Te (x):
    
    psi = Psi(x)

    return mtanh (psi, Teped, Tesep, Tesol, Wd)

# .......................
# Ion temperature profile
# .......................
def Ti (x):
    
    psi = Psi(x)

    return mtanh (psi, Tiped, Tisep, Tisol, Wd)

# ...............................
# Electron number density profile
# ...............................
def ne (x):
    
    psi = Psi(x)

    return mtanh (psi, neped, nesep, nesol, Wd)

# .............    
# tau_R profile
# .............
def tauR (x):

    n_e = ne(x)
    T_e = Te(x) * sc.e
    
    Lambda = 24. - math.log ((n_e/1.e6)**0.5 / (T_e /sc.e))
    tauei  = 6.*2.**0.5*math.pi**1.5 * sc.epsilon_0**2. * sc.m_e**0.5 * T_e**1.5 /Z /Lambda /sc.e**4 /n_e
    etap   = sc.m_e /1.96 /n_e /sc.e**2 /tauei

    return sc.mu_0 * a*a * r*r /etap

# .............
# tau_A profile
# .............
def tauA (x):

    n_e = ne(x)

    return (R0 /B0) * (sc.mu_0 * M * sc.m_p * n_e)**0.5

# ................
# tau_perp profile
# ................
def tauperp (x):

    return a*a * r*r /chiperp

# ...............
# tau_phi profile
# ...............
def tauphi (x):

    return a*a * r*r /chiphi

# ..................
# hat_d_beta profile
# ..................
def dbeta (x):

    T_e = Te(x) * sc.e
    T_i = Ti(x) * sc.e

    return ((5./3.) * M * sc.m_p * (T_e + T_i))**0.5 /sc.e /B0/ a /r

# ...................
# omega_ast_e profile
# ...................
def omegaste (x):

    psip = B0 * a * a * math.exp (q0 /alpham) * sp.exp1 (q0 /alpham) /2. /alpham

    psi = Psi(x)
    T_e = Te(x) * sc.e
    n_e = ne(x)

    dTedp = mtanhp (psi, Teped, Tesep, Tesol, Wd) * sc.e
    dnedp = mtanhp (psi, neped, nesep, nesol, Wd)

    return (dTedp + dnedp * T_e /n_e) /sc.e /psip

# ...................
# omega_ast_i profile
# ...................
def omegasti (x):

    psip = B0 * a * a * math.exp (q0 /alpham) * sp.exp1 (q0 /alpham) /2. /alpham

    psi = Psi(x)
    T_i = Ti(x) * sc.e
    n_e = ne(x)

    dTidp = mtanhp (psi, Tiped, Tisep, Tisol, Wd) * sc.e
    dnedp = mtanhp (psi, neped, nesep, nesol, Wd)

    return - (dTidp + dnedp * T_i /n_e) /sc.e /psip

# ...............
# omega_E profile
# ...............
def omegaE (x):

    we = omegaste (x)

    return omegaEfac * we

# .........
# S profile
# .........
def S (x):

    tau_R = tauR(x)
    tau_A = tauA(x)

    return tau_R /tau_A

# ...........
# Q_e profile
# ...........
def Qe (x):

    s     = S(x)
    tau_A = tauA(x)
    we    = omegaste(x)
    
    return - pow (s, 1./3.) * ntor * we * tau_A

# ...........
# Q_i profile
# ...........
def Qi (x):

    s     = S(x)
    tau_A = tauA(x)
    wi    = omegasti(x)
    
    return - pow (s, 1./3.) * ntor * wi * tau_A

# ...........
# Q_E profile
# ...........
def QE (x):

    s     = S(x)
    tau_A = tauA(x)
    wE    = omegaE(x)
    
    return - pow (s, 1./3.) * ntor * wE * tau_A

# .........
# D profile
# .........
def D(x):

    s   = S(x)
    tau = - omegaste(x) /(omegasti(x) + 1.e-15)
    db  = dbeta(x)

    return pow (s, 1./3.) * (tau /(1. + tau))**0.5 * db

# .............
# Pperp profile
# .............
def Pperp(x):

    tr = tauR(x)
    tp = tauperp(x)

    return tr /tp

# .............
# Pphi profile
# .............
def Pphi(x):

    tr = tauR(x)
    tp = tauphi(x)

    return tr /tp

# ...............
# epsilon profile
# ...............
def epsilon(x):

    if x < 1.:
        return (1. - x*x) /2./ntor /x
    else:
        return (x*x - 1.) /2./ntor /x

# ...............
# Delta profile
# ...............
def Delta(x):

    wE = omegaE(x)
    we = omegaste(x)

    taua = tauA(x)
    taur = tauR(x)
    taup = tauperp(x)
    db   = dbeta(x)
    ns   = ntor * abs (shear(x))

    return (2.*math.pi*math.gamma(0.75) /math.gamma(0.25)) * ntor * abs (wE + we) * taua**0.5 * taur**0.75 /taup**0.25 /db**0.5 /ns**0.5

# ...............
# delta profile
# ...............
def delta(x):

    taua = tauA(x)
    taur = tauR(x)
    taup = tauperp(x)
    db   = dbeta(x)
    ns   = ntor * abs (shear(x))

    return taua**0.5 / taur**0.25 /taup**0.25 /db**0.5 /ns**0.5

# ...........
# Ekk profile
# ...........
def Ekk(x):

    qval = q(x)

    return 2. * ntor * qval

# ################
# Generate figures
# ################

rrr = np.linspace (rmin, rmax, 1000)

ee = []
d1 = []
d2 = []
d3 = []
ps = []
db = []
qx = []

for r in rrr:    

    if r == 1.:
        r = 1. + 1.e-15
        
    ee.append (math.log10(epsilon(r)))
    d1.append (math.log10(delta(r)))
    d2.append (Delta(r))
    d3.append (Delta(r)/Ekk(r))
    ps.append (Psi(r) - 1)
    db.append (math.log10(dbeta(r)))
    qx.append (q(r))

for n in range(len(qx) - 1):

    if (d1[n] - ee[n]) * (d1[n+1] - ee[n+1]) < 0.:
        print ("Psi = %11.4e q = %11.4e m = %3d" % (1.+ps[n], qx[n], int(ntor*qx[n]))) 

fig = plt.figure (figsize = (5.0, 7.0))

plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (3, 1, 1)

plt.xlim (ps[0], ps[-1])
 
plt.plot (ps, db, color = 'green', linewidth = 2, linestyle = 'dashed',  label = r'$\log_{10}(\hat{d}_\beta)$')
plt.plot (ps, ee, color = 'red',   linewidth = 2, linestyle = 'solid',   label = '$log_{10}(\hat{\epsilon}_k)$')
plt.plot (ps, d1, color = 'blue',  linewidth = 2, linestyle = 'solid',   label = '$\log_{10}(\hat{\delta}_k)$')

plt.axvline (0., color = 'black', linewidth = 2, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'x', scilimits = (0, 0))

plt.xlabel (r'$\Psi-1$', fontsize = "15")

plt.legend (fontsize = '12')

plt.subplot (3, 1, 2)

plt.xlim (ps[0], ps[-1])

plt.ylim (0., 30.)

plt.plot (ps, d2, color = 'green',  linewidth = 2, linestyle = 'solid')

plt.axvline (0., color = 'black', linewidth = 2, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'x', scilimits = (0, 0))

plt.xlabel (r'$\Psi-1$',   fontsize = "15")
plt.ylabel (r'$|\Delta_k|$', fontsize = "15")

plt.subplot (3, 1, 3)

plt.xlim (ps[0], ps[-1])

plt.ylim (0., 2.5)

plt.plot (ps, d3, color = 'red',  linewidth = 2, linestyle = 'solid')

plt.axvline (0., color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axhline (1., color = 'black', linewidth = 1, linestyle = 'dotted')

plt.ticklabel_format (style = 'sci', axis = 'x', scilimits = (0, 0))

plt.xlabel (r'$\Psi-1$',             fontsize = "15")
plt.ylabel (r'$|\Delta_k|/(-E_{kk})$', fontsize = "15")
      
plt.tight_layout ()

plt.show ()
#plt.savefig ("Figure9.pdf")
