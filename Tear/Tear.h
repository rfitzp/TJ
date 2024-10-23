// Tear.h

// #########################################################################################

// Class to solve cylindrical tearing mode problem

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

// Inputs:
//  Inputs/Tear.json - JSON file

// Outputs:
//  Outputs/Tear/Tear.nc

// #########################################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>  

#include <nlohmann/json.hpp>
#include <netcdf>

using namespace std;
using           json = nlohmann::json;
using namespace netCDF;
using namespace netCDF::exceptions;
    
// ############
// Class header
// ############
class Tear
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  int    MPOL;   // Poloidal mode number (read from JSON file)
  int    NTOR;   // Toroidal mode number (read from JSON file)
  double mpol;   // Poloidal mode number
  double ntor;   // Toroidal mode number
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
 
  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;     // Integration accuracy (read from JSON file)
  double h0;      // Initial integration step-length (read from JSON file)
  double hmin;    // Minimum integration step-length (read from JSON file)
  double hmax;    // Maximum integration step-length (read from JSON file)
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Integration error calculation flag

  // ----------------
  // Calculation data
  // ----------------
  double qs;     // Resonant safety-factor value
  double rs;     // Radius of rational surface
  double Delta;  // Tearing stability index

  double* rr;     // Radial grid
  double* qq;     // Safety factor
  double* ss;     // Magnetic shear
  double* JJ;     // Toroidal plasma current
  double* JJp;    // Toroidal plasma current gradient
  double* lvals;  // Tearing mode drive term
  double* Psi;    // Tearing mode eigenfunction

  // ----------------------------
  // Cash-Karp RK4/RK5 parameters
  // ----------------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int count;

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

  // Rind rational surface radius
  double FindRationalSurface ();

  // Calculate tearing stability index
  double GetDelta ();
  
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, double* y, double* dydx);

  // Adaptive step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag, 
			   int diag, FILE* file);
  // Fixed step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h);

  // Read JSON file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);
};
