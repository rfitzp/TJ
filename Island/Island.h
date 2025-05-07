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

  // .......................
  // Island calculation data
  // .......................
  double          X;       // Normalized radial coordinate
  double*         XX;      // Normalized radial grid points
  double*         kk;      // k grid points
  double*         zz;      // Angular grid points
  double*         dTdk;    // dT/dk function
  double*         F;       // Temperature perturbation flux-function
  double          T0pls;   // Asymptotic value of X - deltaT[0] at large positive X
  double          T0min;   // Asymptotic value of X + deltaT[0] at large negative X
  double          T0inf;   // Sum of T0pls and T0min
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

  // ...............
  // ECCD parameters
  // ...............
  int     ECCD;   // Flag for eccd calculation (read from JSON file) (0/+1/-1 no calculation/W scan/D scan
  int     Nk;     // Number of k grid points (read from JSON file)
  int     Nscan;  // Number of points in W scan (read from JSON file)
  double  Kmax;   // Maximum value of k on k grid (read from JSON file)
  double  Wmax;   // Maximum value of W in W scan (read from JSON file)
  double  Dmax;   // Maximum value of D in W scan (read from JSON file)
  double  D;      // Radial offset relative to radial width of eccd deposition in W scan (read from JSON file)
  double  W;      // Island width relative to radial width of eccd deposition in D scan (read from JSON file)

  // .....................
  // ECCD calculation data
  // .....................
  double*            kkk;         // k grid points
  double*            Flux0;       // Flux-surface average of 1
  double*            Flux1;       // Flux-surface average of cos(zeta)
  double*            Flux2;       // k <cos(zeta)> /<1>
  double*            Flux3;       // Flux surface average of cos(xi)
  double*            Flux4;       // Flux surface average of sin(xi)*sin(zeta)
  double*            Flux5;       // Flux surface average of Y^2
  double*            Flux6;       // k (<cos xi> + delta^2 <sin xi sin zeta>) <cos zeta> /<1>
  double*            Flux7;       // k <cos zeta> /<1> /<Y^2>
  double*            WW;          // Island width grid
  double*            DD;          // Radial offset grid
  Array<double,2>    JO;          // Flux-surface average of eccd when aimed at O-point
  Array<double,2>    JX;          // Flux-surface average of eccd when aimed at X-point
  Array<double,2>    IO;          // Integrand for Delta_eccd calculation when eccd aimed at O-point
  Array<double,2>    IX;          // Integrand for Delta_eccd calculation when eccd aimed at X-point
  double*            G1;          // Integrand for Delta_ruth calculation
  double*            G2;          // Integrand for Delta_boot calculation

  gsl_spline**       IO_spline;   // Interpolated IO functions
  gsl_spline**       IX_spline;   // Interpolated IX functions
  gsl_spline*        G1_spline;   // Interpolated G1 function
  gsl_spline*        G2_spline;   // Interpolated G2 function
 
  gsl_interp_accel** IO_acc;      // Accelerator for interpolated IO functions
  gsl_interp_accel** IX_acc;      // Accelerator for interpolated IX functions
  gsl_interp_accel*  G1_acc;      // Accelerator for interpolated G1 function
  gsl_interp_accel*  G2_acc;      // Accelerator for interpolated G2 function

  double  DeltaR;                 // Delta_ruth value
  double  DeltaB;                 // Delta_boot value
  double* DeltaO;                 // Delta_eccd values when eccd aimed at O-point
  double* DeltaX;                 // Delta_eccd values when eccd aimed at X-point

  // Misc
  char   buffer[100];     // Name of netcdf file
  double kval, Wval;
  int    rhs_chooser;

public:

  // Constructors
  Island ();
  Island (int k, double _delta);

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

  // Get JO (s, k, xi, W, D)
  double GetJO (int s, double k, double xi, double Wval, double Dval);
  // Get JO_plus (k, xi, W, D)
  double GetJOplus (double k, double xi, double Wval, double Dval);
  // Get JX (s, k, xi, W, D)
  double GetJX (int s, double k, double xi, double Wval, double Dval);
  // Get JX_plus (k, xi, W, D)
  double GetJXplus (double k, double xi, double Wval, double Dval);
  
  // Write Island data to netcdf file
  void WriteNetcdf (int FLAG);
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;

  // Target function for one-dimensional root finding
  double RootFindF (double xi) override;
};
