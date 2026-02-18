// Kink.h

// #########################################################################################

// Class to solve for external kink stabiltiy in cylindrical plasma.

// All lengths normalized to a (minor radius of plasma).
// So r = 0 is magnetic axis and r = 1 is plasma/vacuum interface.
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).

// Equilibrium profiles:

//  Current profile is
//
//  sigma(r) = [2 (nu+1) /qa] (1 - r^2)^nu
//
// qa is safety-factor at plasma/vacuum interface. 
// nu is current peakedness parameter.

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Inputs:
//  Inputs/Kink.json - JSON file

// Outputs:
//  Outputs/Kink/Kink.nc

// #########################################################################################

#pragma once

#include "Utility.h"
    
// ############
// Class header
// ############
class Kink : private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  int    NTOR;   // Toroidal mode number (read from JSON file)
  int    MPOL;   // Poloidal mode number (read from JSON file)
  double qa;     // Edge safety-factor (read from JSON file)
  double nu;     // Current peaking factor (read from JSON file)

  // ---------------
  // Scan parameters
  // ---------------
  int Nqa;    // Number of points in qa scan (read from JSON file)
  int Nnu;    // Number of points in nu scan (read from JSON file)

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from JSON file)
  double del;    // Distance of closest approach to rational surface (read from JSON file)
  int    Nr;     // Number of grid-points (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double qs;     // Resonant safety-factor value
  double rs;     // Radius of rational surface
  double Lambda; // Kink stability parameter
  double bcrit;  // Critical wall radius

  double mpol;   // Poloidal mode number
  double ntor;   // Toroidal mode number
  
  double* rr;    // Radial grid
  double* qq;    // Safety factor
  double* ss;    // Sigma profile
  double* ppsi;  // Kink eigenfunction

  double* qqa;   // Edge safety-factor
  double* ll;    // Lambda
  double* gg;    // gamma
  double* bb;    // b_crit
  double* ww;    // delta W_nw

  double* nn;    // Current peakedness parameter
  double* qxa;   // Edge safety-factor

public:

  // Constructor
  Kink ();
  // Destructor
  ~Kink ();

  // Solve problem
  void Solve ();

private:

  // Find critical qa for marginal stability
  void FindCritical ();
  
  // Scan qa
  void Scanqa ();
  
  // Write data to netcdf file
  void WriteNetcdf ();

  // Return equilibrium quantities
  void GetEquilibrium (double r, double& q, double& sigma, double& drive);

  // Return value of safety-factor
  double Getq (double r);

  // Return value of derivative of safety-factor
  double Getqp (double r);

  // Find rational surface radius
  double FindRationalSurface ();

  // Calculate kink stability index
  double GetLambda ();
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
  // Target function for zero finding
  double RootFindF (double x) override;
};
