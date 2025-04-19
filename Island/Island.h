// Island.h

// #####################################################################

// Class to calculate temperature perturbation in vicinity of asymmetric
// magnetic island

// Inputs:
//  Inputs/Island.json - Island JSON file
//  Inputs/Layer.json  - JSON file

// Outputs:
//  Outputs/Island/Island.nc

// Plotting scripts:
//  Plots/island/*.py

// Class uses following external libraries:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Author:
// Richard Fitzpatrick,
// Institute of Fusion Studies,
// Department of Physics
// University of Texas at Austin
// rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/Island.pdf

// #####################################################################

#pragma once

#include <blitz/array.h>
#include <netcdf>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

#include "Utility.h"

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

// ############
// Class header
// ############
class Island : private Utility
{
private:

  // .................
  // Island parameters
  // .................
  int    Nh;    // Number of harmonics in calculation (read from JSON file)
  int    Nb;    // Number of Bessel harmonics in calculation (read from JSON file)
  int    Nz;    // Number of angular grid points (read from JSON file)
  int    NX;    // Number of radial grid points (read from JSON file)
  double Xmax;  // Maximum value of radial variable (in island widths) (read from JSON file)
  double delta; // Asymmetry parameter (read from JSON file)

  // ................
  // Calculation data
  // ................
  double          X;       // Normalized radial coordinate
  double*         XX;      // Normalized radial grid points
  double*         kk;      // k grid points
  double*         zz;      // Angular grid points
  double*         dTdk;    // dT/dk function
  double*         F;       // Temperature perturbation flux-function
  double          T0inf;   // Asymptotic value of X - deltaT[0] at large X
  double*         dTo;     // Temperature perturbation versus X at zeta = pi
  double*         dTx;     // Temperature perturbation versus X at zeta = 0
  double*         ximx;    // Maximum value of xi

  Array<double,2> En;      // En(1/k) functions
  Array<double,2> deltaTh; // Harmonics of temperature perturbation versus X
  Array<double,2> deltaT;  // Temperature perturbation versus X and zeta

  gsl_spline**        Enspline;   // Interpolated E_n(k) functions
  gsl_spline*         dTdkspline; // Interpolated dTdk(k) function    
  gsl_spline*         Fspline;    // Interpolated F(k) function
  gsl_interp_accel**  Enacc;      // Accelerator for interpolated E_n(k) functions
  gsl_interp_accel*   dTdkacc;    // Accelerator for interpolated dTdk(k) function
  gsl_interp_accel*   Facc;       // Accelerator for interpolated F(k) function

  // Misc
  int rhs_chooser;

public:

  // Constructors
  Island ();
  Island (double _delta);

  // Solve problem
  void Solve (int FLAG);

private:

  // Get G(k)
  double GetG (double k);
  // Get cos(zeta)
  double GetCosZeta (double xi);
  // Get sin(zeta)
  double GetSinZeta (double xi);
  // Get cos(n*zeta)
  double GetCosnZeta (int n, double xi);
  // Get sin(n*zeta)
  double GetSinnZeta (int n, double xi);
  // Get xi(zeta)
  double GetXi (double zeta);
  // Get zeta(xi)
  double GetZeta (double xi);
  // Get xic(X)
  double GetXic (double x);
  // Get k(X, zeta)
  double Getk (double x, double zeta);
  // Get sigma(xi)
  double GetSigma (double xi);
  // Get kappa(x, xi)
  double GetKappa (double x, double xi);
  
  // Write Island data to netcdf file
  void WriteNetcdf (int FLAG);

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;
};
