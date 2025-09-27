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

#include "Utility.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

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
  double qa;      // Safety-factor at plasma boundary
  double beta0;   // Central plasma beta (read from JSON file)
  double alphas;  // Parallel current profile parameter (read from JSON file)
  double nus;     // Parallel current profile parameter (read from JSON file)
  double alphap;  // Pressure profile parameter (read from JSON file)
  double nup;     // Pressure profile parameter (read from JSON file)

  double bwall;   // Minor radius of resistive wall relative to that of plasma (read from JSON file)
  double dwall;   // Radial thickness of resistive wall relative to plasma minor radius (read from JSON file)
  double epsb;    // Inverse aspect-ratio of resistive wall
  double delw;    // Wall thickness parameter

  double gmax;    // Maximum growth/decay rate (read from JSON file)
  
  double Theta;   // Pinch parameter 
  double Frev;    // Reversal parameter 

  // ---------------
  // Mode parameters
  // ---------------
  int    mpol;   // Poloidal mode number (read from JSON file)
  int    ntor;   // Toroidal mode number (read from JSON file)
  double MPOL;   // Poloidal mode number
  double NTOR;   // Toroidal mode number
  double ka;     // Toroidal wavenumber at plasma boundary
  double kb;     // Toroidal wavenumber at wall
  double qres;   // Resonant safety-factor
  double rres;   // Radius of resonant surface
  
  // ----------------------
  // Calculation parameters
  // ----------------------
  int    Ngrid;  // Number of radial grid points (read from JSON file)
  double eps;    // Closest approach to magnetic axis (read from JSON file)
  double delta;  // Closest approach to rational surface (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double* rr;       // Radial grid
  double* rrv;      // Radial grid in vacuum
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

  double* bbeta1;   // r P' profile
  double* bbeta2;   // r^2 P'' profile
  double* ssigp;    // sigma' profile
  double* ff;       // f profile
  double* gg;       // tanh(g/10) profile

  gsl_interp_accel* Bphi_accel;    // Interpolation accelerator for toroidal magnetic field
  gsl_spline*       Bphi_spline;   // Interpolator for toroidal magnetic field
  gsl_interp_accel* Btheta_accel;  // Interpolation accelerator for poloidal magnetic field
  gsl_spline*       Btheta_spline; // Interpolator for poloidal magnetic field
  gsl_interp_accel* Pp_accel;      // Interpolation accelerator for pressure gradient
  gsl_spline*       Pp_spline;     // Interpolator for pressure gradient
  gsl_interp_accel* q_accel;       // Interpolation accelerator for safety-factor
  gsl_spline*       q_spline;      // Interpolator for safety-factor

  double* Psip;    // Plasma solution
  double* Psinw;   // No-wall vacuum solution
  double* Psipw;   // Perfect-wall vacuum solution
  double* Psirwm;  // Resistive wall mode vacuum solution

  double lbar;  // Plasma stability index
  double Lnw;   // No-wall vacuum stability index
  double Lpw;   // No-wall vacuum stability index
  double alpw;  // Wall parameter
  double Wnw;   // No-wall delta-W
  double Wpw;   // Perfect-wall delta-W
  double c1;    // Amount of Psinw in Psirwm
  double c2;    // Amount of Psipw in Psirwm
  double rhs;   // Right-hand side of resistive wall mode dispersion relation

  double gamma; // Resistive wall mode growth-rate

  // ----
  // Misc
  // ----
  int count, rhs_chooser, f_chooser;
  
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
  // Find resonant surface
  void FindResonant ();
  // Solve Newcomb's equation
  void SolveNewcomb ();
  // Solve vacuum solution
  void SolveVacuum ();
  // Calculate resistive wall mode growth-rate
  void Growth ();
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
  // Get q
  double Getq (double r);
  // Get Btheta
  double GetBtheta (double r);
  // Get Bphi
  double GetBphi (double r);
  // Get beta1
  double Getbeta1 (double r);
  // Get beta2
  double Getbeta2 (double r);
  // Get F
  double GetF (double r);
  // Get G
  double GetG (double r);
  // Get H
  double GetH (double r);
  // Get f
  double Getf (double r);
  // Get g
  double Getg (double r);
  // Get Im
  double GetIm (double r);
  // Get Km
  double GetKm (double r);
  // Get Imp
  double GetImp (double r);
  // Get Km
  double GetKmp (double r);
 
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;
  // Target function for zero finding
  double RootFindF (double x) override;
};

