// Island.h

// #####################################################################

// Class to calculate temperature perturbation in vicinity of magnetic
// island

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
  int    Nz;    // Number of angular grid points (read from JSON file)
  int    Nx;    // Number of radial grid points (read from JSON file)
  double xmax;  // Maximum value of radial variable (in island widths) (read from JSON file)

  // ................
  // Calculation data
  // ................
  double          x;       // Radial coordinate
  double*         xx;      // Radial grid points
  double*         kk;      // k grid points
  double*         zz;      // Angular grid points
  double*         F;       // Temperature perturbation flux-function
  double          Finf;    // Asymptotic value of x - F at large x
  double*         dTo;     // Temperature perturbation versus x at zeta = pi
  double*         dTx;     // Temperature perturbation versus x at zeta = 0
  Array<double,2> deltaTh; // Harmonics of temperature perturbation versus x
  Array<double,2> deltaT;  // Temperature perturbation versus x and zeta
  
  gsl_spline*        Fspline;  // Interpolated F(k) function
  gsl_interp_accel*  Facc;     // Accelerator for interpolated F(k) function

  // Misc
  int rhs_chooser;

public:

  // Constructor
  Island ();

  // Solve problem
  void Solve ();

private:

  // Get k
  double Getk (double x, double zeta);
  // Get zeta_c
  double Getzetac (double x);
  
  // Write Island data to netcdf file
  void WriteNetcdf ();

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;
};
