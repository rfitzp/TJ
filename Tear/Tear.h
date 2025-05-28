// Tear.h

// #########################################################################################

// Class to solve cylindrical tearing mode problem. Class calculates Delta' for all rational
// surfaces in plasma.

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

// Inputs:
//  Inputs/Tear.json - JSON file

// Outputs:
//  Outputs/Tear/Tear.nc

// #########################################################################################

#pragma once

#include <blitz/array.h>
#include <netcdf>
#include "Utility.h"

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
    
// ############
// Class header
// ############
class Tear : private Utility
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

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from JSON file)
  double del;    // Distance of closest approach to rational surface (read from JSON file)
  int    Nr;     // Number of grid-points (read from JSON file)

   // ----------------
  // Calculation data
  // ----------------
  double          qs;    // Resonant safety-factor value
  double          rs;    // Radius of rational surface
  double          Delta; // Tearing stability index

  double          mpol;  // Poloidal mode number
  double          ntor;  // Toroidal mode number
  int             nres;  // Number of rational surfaces
  int*            mres;  // Resonant poloidal mode numbers
  double*         rres;  // Rational surface radii
  double*         Dres;  // Tearing stability indicies
  
  double*         rr;    // Radial grid
  double*         qq;    // Safety factor
  double*         ss;    // Magnetic shear
  double*         JJ;    // Toroidal plasma current
  double*         JJp;   // Toroidal plasma current gradient
  double*         lvals; // Tearing mode drive term
  Array<double,2> Psi;   // Tearing mode eigenfunctions

public:

  // Constructor
  Tear ();
  // Destructor
  ~Tear ();

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
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
};
