// Pinch.h

// ###############################################################
// Class to calculate toroidal pinch equilibria

// Parallel current profile:
//
//  sigma(r) = (2*epsa/q0) (1 - r^alphas)^nus
//
// Pressure profile:
//
// mu_0 * P(r) /[B_phi(0)]^2 = (beta0/2) (1 - r^alphap)^nup
//

// Inputs:
//  Inputs/Pinch.json - JSON file

// Outputs:
//  Outputs/Pinch/Pinch.nc

// Plotting scripts:
//  Plots/Pinch/*.py

// Class uses following external libraries:
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics
//  University of Texas at Austin
//  rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ
// ###############################################################

#pragma once

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <netcdf>
#include "Utility.h"

using namespace netCDF;
using namespace netCDF::exceptions;
    
// ############
// Class header
// ############
class Pinch : private Utility
{
 private:

  // ----------------------
  // Equilibrium parameters
  // ----------------------
  double epsa;    // Inverse aspect-ratio of plasma boundary (read from JSON file)
  double q0;      // Central safety-factor (read from JSON file)
  double beta0;   // Central plasma beta (read from JSON file)
  double alphas;  // Parallel current profile parameter (read from JSON file)
  double nus;     // Parallel current profile parameter (read from JSON file)
  double alphap;  // Pressure profile parameter (read from JSON file)
  double nup;     // Pressure profile parameter (read from JSON file)

  double Theta;   // Pinch parameter 
  double Frev;    // Reversal parameter 
  
  // ----------------------
  // Calculation parameters
  // ----------------------
  int    Ngrid;  // Number of radial grid points (read from JSON file)
  double eps;    // Closest approach to magnetic axis (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double* rr;       // Radial grid
  double* ssigma;   // Parallel current profile 
  double* PP;       // Pressure profile 
  double* BBphi;    // Toroidal magnetic field 
  double* BBtheta;  // Poloidal magnetic field 
  double* qq;       // Safety-factor profile

  double* qqc;      // Mercier-stable safety-factor profile
  double* PPc;      // Mercier-stable pressure profile
  double* BBphic;   // Mercier-stable toroidal magnetic field 
  double* BBthetac; // Mercier-stable poloidal magnetic field

  double* PPp;      // Pressure gradient profile
  double* PPpc;     // Critical pressure gradient profile
  double* PPpm;     // Mercier-stable pressure gradient profile

  // ----
  // Misc
  // ----
  int count, rhs_chooser;
  
public:

  // Constructor
  Pinch ();
  // Destructor
  ~Pinch ();

  // Set beta0
  void Setbeta0 (double _beta0);
  // Set q0
  void Setq0 (double _q0);
  // Set nus
  void Setnus (double _nus);
  // Solve problem
  void Solve ();

private:

  // Calculate equilibrium
  void CalcEquilibrium ();
  // Output data to Netcdf file
  void WriteNetcdf ();
  
  // Get equilibrium parallel current
  double GetSigma (double r);
  // Get equilibrium parallel current gradient
  double GetSigmap (double r);
  // Get equilibrium pressure
  double GetP (double r);
  // Get equilibrium pressure gradient
  double GetPp (double r);
  // Get magnetic shear
  double Gets (double r, double q);
  // Get critical pressure gradient
  double GetPpcrit (double r, double q, double Bphi);

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;
};

