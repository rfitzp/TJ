// TearX.h

// #########################################################################################

// Class to solve cylindrical tearing mode problem for equilibrium with simulated X-point.
// Class calculates Delta' for all rational surfaces in plasma associated with given
// toroidal mode number.

// All lengths normalized to a (minor radius of plasma).
// So r = 0 is magnetic axis and r = 1 is plasma/vacuum interface.
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).

// Equilibrium profiles:

//  Safety-factor profile is q(r) = r^2 /f(r)
//
//  f(r) = (1 /nu/q0) [1 - (1 - r^2)^nu] - alpha ln(1 - r^2)
//
// q0 is safety-factor on magnetic axis.
// qa = nu * q0 is cylindrical safety-factor at plasma/vacuum interface.

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
    
// ############
// Class header
// ############
class TearX: private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  int    NTOR;   // Toroidal mode number (read from JSON file)
  double q0;     // Central safety-factor (read from JSON file)
  double nu;     // Current peaking factor (read from JSON file)
  double qc;     // Cylindrical edge safety-factor 
  double alpha;  // X-point parameter (read from JSON file)

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from JSON file)
  double del;    // Distance of closest approach to rational surface (read from JSON file)
  int    Nr;     // Number of radial grid-points (read from JSON file)
  double EPS;    // ln(1 - r^2) regularized as ln(1 - r^2 + EPS) (read from JSON file)
  double Psimax; // Rational surfaces ignored in region Psi > Psimax (read from JSON file)

  double q0_sta; // Start value for q0 scan  (read from JSON file)
  double q0_end; // End value for q0 scan  (read from JSON file)
  int    q0_num; // Number of points in q0 scan  (read from JSON file)
  double nu_sta; // Start value for nu scan  (read from JSON file)
  double nu_end; // End value for nu scan  (read from JSON file)
  int    nu_num; // Number of points in nu scan  (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double            r95;      // Radius of 95% flux-surface
  double            q95;      // Safety-factor at 95% flux-surface
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
  double*           sres;     // Magnetic shears
  double*           Dres;     // Tearing stability indicies
  
  double*           rr;       // Radial grid
  double*           qq;       // Safety factor
  double*           PSI;      // Normalized equilibrium poloidal magnetic flux
  double*           ss;       // Magnetic shear
  double*           JJ;       // Toroidal plasma current
  double*           JJp;      // Toroidal plasma current gradient
  double*           lvals;    // Tearing mode drive term
  Array<double,2>   Psi;      // Tearing mode eigenfunctions

  gsl_spline*       q_spline; // Interpolated q function
  gsl_spline*       P_spline; // Interpolated PSI function

  gsl_interp_accel* q_acc;    // Accelerator for interpolated q function
  gsl_interp_accel* P_acc;    // Accelerator for interpolated PSI function

  // ----
  // Misc
  // ----
  int rhs_chooser;

public:

  // Constructor
  TearX ();
  // Destructor
  ~TearX ();

  // Scan q0
  void Scanq0 ();
  // Scan nu
  void Scannu ();
  
  // Solve problem
  void Solve (int verbose);

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
  // Return value of magnetic shear
  double GetShear (double r);

  // Find rational surface radius
  double FindRationalSurface ();

  // Calculate tearing stability index
  double GetDelta (int isurf);
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
};
