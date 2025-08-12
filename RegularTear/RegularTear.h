// Tear.h

// #########################################################################################

// Class to solve regularized cylindrical tearing mode problem. Class calculates Delta'
// for all rational surfaces in plasma.

// All lengths normalized to a (minor radius of plasma).
// So r = 0 is magnetic axis and r = 1 is plasma/vacuum interface.
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).

// Equilibrium profiles:

//  Safety-factor profile is q(r) = r^2 /f(r)
//
//  f(r) = (1 /nu/q0) [1 - (1 - r^2)^nu] 
//
// q0 is safety-factor on magnetic axis.
// qa = nu * q0 is safety-factor at plasma/vacuum interface.

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)
//  gsl library            (https://https://www.gnu.org/software/gsl)

// Inputs:
//  Inputs/RegularTear.json - JSON file

// Outputs:
//  Outputs/RegularTear/RegularTear.nc

// #########################################################################################

#pragma once

#include <gsl/gsl_sf_erf.h>

#include "Utility.h"
    
// ############
// Class header
// ############
class RegularTear : private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  int    NTOR;   // Toroidal mode number (read from JSON file)
  double q0;     // Central safety-factor (read from JSON file)
  double qa;     // Edge safety-factor (read from JSON file)
  double nu;     // Current peaking factor
  int    Fixed;  // Flag for fixed boundary calculation
  double sigma;  // Regularization parameter

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from JSON file)
  double del;    // Distance of closest approach to rational surface (read from JSON file)
  int    Nr;     // Number of grid-points (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double          qs;     // Resonant safety-factor value
  double          rs;     // Radius of rational surface
  double          sh;     // Magnetic shear at rational surface
  double          Delta;  // Tearing stability index
  double          Deltr;  // Regularized tearing stability index

  double          mpol;   // Poloidal mode number
  double          ntor;   // Toroidal mode number
  int             nres;   // Number of rational surfaces
  int*            mres;   // Resonant poloidal mode numbers
  double*         rres;   // Rational surface radii
  double*         sres;   // Magnetic shear at rational surfaces
  double*         Dres;   // Tearing stability indicies
  double*         Dresr;  // Regularized tearing stability indicies
  
  double*         rr;     // Radial grid
  double*         qq;     // Safety factor
  double*         ss;     // Magnetic shear
  double*         JJ;     // Toroidal plasma current
  double*         JJp;    // Toroidal plasma current gradient
  double*         lvals;  // Tearing mode drive term

  Array<double,2> Psi;    // Tearing mode eigenfunctions
  Array<double,2> Psip;   // Derivatives of tearing mode eigenfunctions
  Array<double,2> Psipp;  // Second derivatives of tearing mode eigenfunctions
  Array<double,2> Cur;    // Tearing mode current eigenfunctions
  Array<double,1> Norm;   // Normalizing factors inside rational surfaces
  Array<double,1> Norp;   // Normalizing factors outside rational surfaces

  Array<double,2> Psir;   // Regularized tearing mode eigenfuctions
  Array<double,2> Psipr;  // Regularized derivatives of tearing mode eigenfunctions
  Array<double,2> Psippr; // Regularized second derivatives of tearing mode eigenfunctions
  Array<double,2> Curr;   // Regularized tearing mode current eigenfunctions

public:

  // Constructor
  RegularTear ();
  // Destructor
  ~RegularTear ();

  // Solve problem
  void Solve ();

private:

  // Write data to netcdf file
  void WriteNetcdf ();

  // Return equilibrium quantities
  void GetEquilibrium (double r, double& q, double& s, double& J, double& Jp, double& lambda);

  // Return value of safety-factor
  double Getq (double r);

  // Return value of derivative of safety-factor
  double Getqp (double r);

  // Find rational surface radius
  double FindRationalSurface ();

  // Calculate tearing stability index
  double GetDelta (int isurf);

   // Calculate regularized tearing stability index
  double GetRegularDelta (int isurf);

  // Return fit for Psi
  double GetPsiFit (double X, double C0 , double C1, double C2, double C3, double C4, double C5);

  // Return fit for Psip
  double GetPsipFit (double X, double C0, double C1, double C2, double C3, double C4);

  // Return fit for Psipp
  double GetPsippFit (double X, double C0, double C1, double C2, double C3);
  
  // Return fit for Current
  double GetCurFit (double X, double C0, double C1, double C2, double C3, double C4, double C5);

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
};
