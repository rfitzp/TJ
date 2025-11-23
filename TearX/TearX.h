// TearX.h

// #########################################################################################

// Class to solve cylindrical tearing mode problem for equilibrium with simulated X-point.
// Class calculates Delta' for all rational surfaces in plasma associated with given
// toroidal mode number.

// Major radius of plasma is R_0.
// All lengths normalized to a (minor radius of plasma).
// So r = 0 is magnetic axis and r = 1 is plasma/vacuum interface.
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).
// Ion mass number is M. Perpendicular energy diffusivity is chi_E. Perpendicular momentum diffusivity is chi_phi.
// Temperatures are in eV. 

// Equilibrium quantities:

//  The safety-factor profile is
//
//   q(r) = r^2 /f(r) - alpha ln(|1 - r^2| + EPS) - (1/2) Delta_q (1 + tanh[(r - r_q) /delta_q]), 
//
//  where
//
//   f(r) = (1 /nu/q0) [1 - (1 - r^2)^nu]).
//
//  The normalized current profile is
//
//   J(r) = (2/q0) (1 - r^2)^(nu-1).
//
//  q0 is safety-factor on magnetic axis.
//  qc = nu * q0 is cylindrical safety-factor at plasma/vacuum interface.
//
//  The equilibrium poloidal magnetic flux satisfies
//
//   dpsi_p/dr = r/q,
//
//  where psi_p(0) = 0. The normalized poloidal flux is PSI(r) = psi_p(r)/psi_p(1).
//
//  The electron number density, electron temperature, and ion temperature profiles are:
//
//   n_e(r) = (ne_0 - ne_1) (1 - r^2)^nu_ne + ne_1 + (1/2) Delta_ne (1 - tanh[(r - r_ne) /delta_ne]),
//   T_e(r) = (Te_0 - Te_1) (1 - r^2)^nu_Te + Te_1 + (1/2) Delta_Te (1 - tanh[(r - r_Te) /delta_Te]),
//   T_i(r) = (Ti_0 - Ti_1) (1 - r^2)^nu_Ti + Ti_1 + (1/2) Delta_Ti (1 - tanh[(r - r_Ti) /delta_Ti]).
//
//   J_bs = (mu_0 R_0^2 e /a^2 B_0^2) (q_cly /r) (r*a/R0)^1/2 (2.40 (T_e + T_i) dn_e/dr + 0.61 n_e dT_e/dr - 0.41 n_e dT_e/dr)
//
//  where q_cyl = r^2 /f(r).
//
//  The electron diamagnetic frequency, ion diamagnetic frequency, and ExB frequency are:
//
//   omega_ast_e =   (T_e dn_e/dr /n_e + dT_e/dr) (q /(B_0 r a^2)),
//   omega_ast_i = - (T_i dn_i/dr /n_i + dT_i/dr) (q /(B_0 r a^2)),
//   omega_E     = alphaE omega_ast_e + omegaE_0 (1 - r^2)^nu_E

// Layer Quantities:
//
//  ln Lambda = 24 - ln [(ne/10^6)^(1/2) / Te]
//  tau_ei    = 6 2^0.5 pi^1.5 epsilon_0^2 m_e^1/2 T_e^1.5 /(Z ln Lambda e^2.5 n_e)
//  fp        = 1. - 1.46 (r*a/R0)^(1/2) + 0.46 * (r*a/R0)^(3/2)
//  eta_par   = m_e /(1.96 ne e^2 tau_ei fp)
//  tau_R     = mu_0 * r*r * a*a /eta_par
//  tau_A     = R_0 * (mu_0 * M * m_p * n_e)^(1/2) / B_0
//  tau_E     = r*r * a*a /chi_E
//  tau_phi   = r*r * a*a /chi_phi
//  beta      = (5./3.) mu_0 n_e (T_e + T_i) e /B0^2
//  c_beta    = (beta /(1 + beta))^(1/2)
//  d_i       = [M m_p /(ne e^2 mu_0)]^(1/2)
//  d_beta    = c_beta d_i
//  S         = tau_R /tau_A
//  iota_e    = - omega_ast_e /(omega_ast_i - omega_ast_e)

//  tau_H     = tau_A /(n |s|)
//  S_H       = S (n |s|)
//  
//  Q_E       = - S^(1/3) n omega_E     tau_A / (n |s|)^(2/3)
//  Q_e       = - S^(1/3) n omega_ast_e tau_A / (n |s|)^(2/3)
//  Q_i       = - S^(1/3) n omega_ast_i tau_A / (n |s|)^(2/3)
//  D         = S^(1/3) (n |s|)^(1/3) iota_e^(1/2) d_beta /a/r
//  P_E       = tau_R /tau_E
//  P_phi     = tau_R /tau_phi
//
//  Delta_s   = Scale Delta
//  Scale     = S^(1/3) (n |s|)^1/3

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Inputs:
//  Inputs/TearX.json - JSON file

// Outputs:
//  Outputs/TearX/TearX.nc

// #########################################################################################

#pragma once

#include "Utility.h"
#include <gsl/gsl_const_mksa.h>
#include "FourField.h"
    
// ############
// Class header
// ############
class TearX: private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  int    NTOR;     // Toroidal mode number (read from JSON file)
  double q0;       // Central safety-factor (read from JSON file)
  double nu;       // Current peaking factor (read from JSON file)
  double alpha;    // X-point parameter (read from JSON file)
  double Delta_q;  // Height of safety-factor step (read from JSON file)
  double r_q;      // Central radius of safery-factor step (read from JSON file)
  double delta_q;  // Width of safety-factor step (read from JSON file)

  double B0;       // Toroidal magnetic field-strength (SI) (read from JSON file)
  double R0;       // Plasma minor radius (SI) (read from JSON file)
  double a;        // Plasma major radius (SI) (read from JSON file)
  double M;        // Ion mass number (read from JSON file)
  double chiE;     // Perpendicular energy diffusivity (SI) (read from JSON file)
  double chip;     // Perpendicular momentum diffusivity (SI) (read from JSON file)
  double Z;        // Effective ion charge number (read from JSON file)

  double ne_0;     // Central electron number density (SI) (read from JSON file)
  double ne_1;     // Electron number density offset (SI) (read from JSON file)
  double Delta_ne; // Height of electron number density pedestal (SI) (read from JSON file)
  double nu_ne;    // Electron number density central peaking parameter (read from JSON file)
  double r_ne;     // Central radius of electron number density pedestal (read from JSON file)
  double delta_ne; // Width of electron number density pedestal (read from JSON file)

  double Te_0;     // Central electron temperature (eV) (read from JSON file)
  double Te_1;     // Rlectron temperature offset (eV) (read from JSON file)
  double Delta_Te; // Height of electron temperature pedestal (eV) (read from JSON file)
  double nu_Te;    // Electron temperature central peaking parameter (read from JSON file)
  double r_Te;     // Central radius of electron temperature pedestal (read from JSON file)
  double delta_Te; // Width of electron temperature pedestal (read from JSON file)

  double Ti_0;     // Central ion temperature (eV) (read from JSON file)
  double Ti_1;     // Ion temperature offset (eV) (read from JSON file)
  double Delta_Ti; // Height of ion temperature pedestal (eV) (read from JSON file)
  double nu_Ti;    // Ion temperature central peaking parameter (read from JSON file)
  double r_Ti;     // Central radius of ion temperature pedestal (read from JSON file)
  double delta_Ti; // Width of ion temperature pedestal (read from JSON file)

  double alphaE;   // ExB frequency factor
  double omegaE_0; // Central ExB frequency offset (rad/s) (read from JSON file)
  double nu_E;     // ExB frequency central peaking parameter (read from JSON file)

  double e;          // Magnitude of electron charge (SI)
  double m_e;        // Mass of electron (SI)
  double m_p;        // Mass of proton (SI)
  double epsilon_0;  // Vacuum permittivity (SI)
  double mu_0;       // Vacuum permeability (SI)

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;      // Distance of closest approach to magnetic axis (read from JSON file)
  double del;      // Distance of closest approach to rational surface (read from JSON file)
  int    Nr;       // Number of radial grid-points (read from JSON file)
  double EPS;      // ln (1 - r^2) regularized as ln (|1 - r^2| + EPS) (read from JSON file)
  double Psimax;   // Rational surfaces ignored in region PSI > Psimax (read from JSON file)

  double nu_sta;   // Start value for nu scan (read from JSON file)
  double nu_end;   // End value for nu scan (read from JSON file)
  int    nu_num;   // Number of points in nu scan (read from JSON file)
  int    m_min;    // Ignore rational surfaces with poloidal mode numbers less than m_min (read from JSON file)
  double r_min;    // Ignore rational surfaces with radii less than r_min (read from JSON file)
  double c_min;    // Use three-field layer Delta calculation for c_beta < c_min (read from JSON file)
  
  // ----------------
  // Calculation data
  // ----------------
  double            qc;       // Cylindrical edge safety-factor 
  double            r95;      // Radius of 95% flux-surface
  double            R95;      // Rho value of 95% flux-surface
  double            q95;      // Safety-factor at 95% flux-surface
  double            s95;      // Magnetic shear (in terms of r) at 95% flux-surface
  double            S95;      // Magnetic shear (in terms of rho) at 95% flux-surface
  double            rmax;     // Radius of PSI = Psimax flux-surface
  double            qmax;     // Safety-factor at PSI = Psimax flux-surface
    
  double            qs;       // Resonant safety-factor value
  double            rs;       // Radius of rational surface
  double            sr;       // Magnetic shear at rational surface
  double            Delta;    // Tearing stability index

  double            mpol;     // Poloidal mode number
  double            ntor;     // Toroidal mode number
  int               nres;     // Number of rational surfaces
  int*              mres;     // Resonant poloidal mode numbers
  double*           rres;     // Rational surface radii
  double*           Pres;     // PSI values
  double*           sres;     // Magnetic shears
  double*           Dres;     // Tearing stability indicies
  double*           Drres;    // Real parts of layer response indices
  double*           Dires;    // Imaginary parts of layer response indices
  double*           Idres;    // Ideal response indicies
  
  double*           rr;       // Radial grid
  double*           qq;       // Safety factor
  double*           PSI;      // Normalized equilibrium poloidal magnetic flux
  double*           rho;      // rho = PSI^1/2
  double*           ss;       // Magnetic shear
  double*           JJ;       // Normalized toroidal plasma current density
  double*           JJp;      // Normalized toroidal plasma current density gradient
  double*           JJbs;     // Normalized bootstrap current density
  double*           lvals;    // Tearing mode drive term
  Array<double,2>   Psi;      // Tearing mode eigenfunctions

  double*           ne;       // Electron number density 
  double*           dnedr;    // Electron number density gradient
  double*           Te;       // Electron temperature
  double*           dTedr;    // Electron temperature gradient
  double*           Ti;       // Ion temperature
  double*           dTidr;    // Ion temperature gradient

  double*           waste;    // Electron diamagnetic frequency
  double*           wasti;    // Ion diamagnetic frequency
  double*           wE;       // ExB frequency

  double*           tauR;     // Resistive diffusion time
  double*           tauA;     // Alfven time
  double*           tauE;     // Energy confinement time
  double*           taup;     // Momentum confinement time
  double*           S;        // Lundquist number
  double*           cbeta;    // Beta parameter
  double*           dbeta;    // Ion sound radius

  double*           QE;       // Normalized ExB frequency
  double*           Qe;       // Normalized electron diamagnetic frequency
  double*           Qi;       // Normalized ion diamagnetic frequency
  double*           D;        // Normalized ion sound radius
  double*           P_E;      // Energy magnetic Prandtl number
  double*           P_phi;    // Momentum magnetic Prandtl number
  double*           Scale;    // Scale factor

  gsl_spline*       q_spline; // Interpolated q function
  gsl_spline*       P_spline; // Interpolated PSI function
  gsl_spline*       R_spline; // Interpolated rho function

  gsl_interp_accel* q_acc;    // Accelerator for interpolated q function
  gsl_interp_accel* P_acc;    // Accelerator for interpolated PSI function
  gsl_interp_accel* R_acc;    // Accelerator for interpolated rho function

  // ----
  // Misc
  // ----
  int rhs_chooser;

public:

  // Constructor
  TearX ();
  // Destructor
  ~TearX ();

  // Scan nu
  void Scannu ();
  
  // Solve problem
  void Solve (int flg);

private:

  // Set q0
  void Setq0 (double _q0);
  // Set nu
  void Setnu (double _nu);

  // Write data to netcdf file
  void WriteNetcdf ();

  // Return equilibrium quantities
  void GetEquilibrium (double r, double& q, double& s, double& J, double& Jp, double& lambda);

  // Return value of safety-factor
  double Getq (double r);
  // Return value of derivative of safety-factor
  double Getqp (double r);
  // Return value of cylindrical safety-factor
  double Getq_cyl (double r);
  // Return value of magnetic shear
  double GetShear (double r);
  // Return value of current density
  double GetJ (double r);
  // Return value of current density gradient
  double GetJp (double r);
  // Return electron number density
  double Getne (double r);
  // Return electron number density gradient
  double Getdnedr (double r);
  // Return electron temperature
  double GetTe (double r);
  // Return electron temperature gradient
  double GetdTedr (double r);
  // Return ion temperature
  double GetTi (double r);
  // Return ion temperature gradient
  double GetdTidr (double r);
  // Return bootstrap current density
  double GetJbs (double r);
  // Return value of electron diamagnetic frequency
  double Getwaste (double r);
  // Return value of ion diamagnetic frequency
  double Getwasti (double r);
  // Return value of ExB frequency
  double GetwE (double r);

  // Return value of Coulomb logarithm
  double GetLambda (double r);
  // Return value of resistive diffusion time
  double GettauR (double r);
  // Return value of Alfven time
  double GettauA (double r);
  // Return value of energy confinement time
  double GettauE (double r);
  // Return value of momentum confinement time
  double Gettaup (double r);
  // Return value of cbeta
  double Getcbeta (double r);
  // Return value of dbeta
  double Getdbeta (double r);

  // Return value of QE
  double GetQE (double r);
  // Return value of Qe
  double GetQe (double r);
  // Return value of Qi
  double GetQi (double r);
  // Return value of D
  double GetD (double r);
  // Return value of P_E
  double GetP_E (double r);
  // Return value of P_phi
  double GetP_phi (double r);
  // Return value of scale factor
  double GetScale (double r);

  // Find rational surface radius
  double FindRationalSurface ();

  // Calculate tearing stability index
  double GetDelta (int isurf);
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
};
