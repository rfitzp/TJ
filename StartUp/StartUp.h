// StartUp.h

// #####################################################################################

// Class to analyze tokamak startup/shutdown

// Inputs:
//  Inputs/StartUp.json - JSON file

// Outputs:
//  Outputs/Startup/StartUp.nc

// Plotting scripts:
//  Plots/StartUp/*.py

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics,
//  University of Texas at Austin,
//  rfitzp@utexas.edu

// Source: https://github.com/rfitzp/StartUp

// ###################################################################################

#pragma once

#include "Utility.h"

// ############
// Class header
// ############
class StartUp : private Utility
{
 private:

  // Physics parameters
  double alpha;  // Parameter that determines normalized perpendicular diffusivity profile (read from JSON file)
                 // chi(x) = f(alpha) (1 + x^2)^alpha
                 // f(alpha) = (1 + alpha) /(2^(1 + alpha) - 1)
  double lambda; // Temperature profile eigenvalue
  double Bt1;    // Bt(1)
  double qa;     // Edge safety-factor profile
  double li;     // Normalized plasma self-inductance
  double Tramp;  // Maximium normalized central electron temperature
  double tramp;  // Minimum normalized ramp time
  double Eramp;  // Normalized electric field
  double Pramp;  // Normalized ohmic heating power

  // Calculation parameters
  double eps;    // Solutions launch from x = eps (read from JSON file)
  double lmin;   // Minimum value lambda in search interval (read from JSON file)
  double lmax;   // Maximum value lambda in search interval (read from JSON file)
  double xmax;   // Maximum value of x (read from JSON file)

  int     Nx;    // Number of radial grid points (read from JSON file)
  double* xx;    // Radial grid points
  double* T;     // Temperature profile
  double* Bt;    // B_theta profile
  double* q;     // Safety-factor profile

  double  amax;   // Maximum value of alpha in alpha scan (read from JSON file)
  int     Na;     // Number of points in alpha scan (read from JSON file)
  double* aa;     // alpha grid points
  double* ll;     // lambda values
  double* qqa;    // qa values
  double* lli;    // li values
  double* TTr;    // Tramp values
  double* ttr;    // tramp values
  double* EEr;    // Eramp values
  double* PPr;    // Pramp values

  // ----
  // Misc
  // ----
  int rhs_chooser;
   
 public:

  // Constructor
  StartUp ();
  // Destructor
  ~StartUp ();

  // Solve problem
  void Solve ();
 
 private:

  // Solve problem for alpha scan
  void SolveScan ();
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;
  
  // Target function for one-dimensional root finding
  double RootFindF (double x) override;

  // Write netcdf file
  void Write_netcdf ();
};
