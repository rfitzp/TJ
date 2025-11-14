// Utility.h

// ############################################################################################

// Class that contains various useful data and functions that can be inherited by other classes

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics
//  University of Texas at Austin
//  rfitzp@utexas.edu

// Class uses following external libraries:
//  nclohmann JSON library             (https://github.com/nlohmann/json)
//  Blitz++ library                    (https://github.com/blitzpp/blitz)
//  GNU scientific library             (https://www.gnu.org/software/gsl)
//  netcdf-c++ library                 (https://github.com/Unidata/netcdf-cxx4)
//  Armadillo library                  (https://arma.sourceforge.net)
//  Steven Johsnson's Faddeeva pacakge (http://ab-initio.mit.edu/faddeeva)

// Source: https://github.com/rfitzp/TJ

// ############################################################################################

#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES
#define ARMA_WARN_LEVEL 0

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef _WIN32
 #include <direct.h>
 #define mkdir _mkdir
#else
 #include <sys/stat.h>
 #include <sys/types.h>
#endif

#include <nlohmann/json.hpp>
#include <blitz/array.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <netcdf>
#include <armadillo>

#include "Faddeeva.h"

using namespace std;
using           json = nlohmann::json;
using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace arma;
//using namespace Faddeeva;

// ############
// Class header
// ############
class Utility
{
private:

  // ----------------------------
  // Cash-Karp RK4/RK5 parameters
  // ----------------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

public:
  
  // -------------------------------------------------
  // Cash-Karp RK4/RK5 adaptive integration parameters
  // -------------------------------------------------
  double acc;     // Integration accuracy 
  double h0;      // Initial integration step-length 
  double hmin;    // Minimum integration step-length 
  double hmax;    // Maximum integration step-length 
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Integration error calculation flag
  int    count;   // Counter for step recalculations 

  // ---------------------------------------
  // One-dimensional root finding parameters
  // ---------------------------------------
  int    Nint;    // Number of search intervals
  double Eta;     // Minimum magnitude of F at root F(x) = 0
  int    Maxiter; // Maximum number of iterations

  // ---------------------------------------
  // Two-dimensional root finding parameters
  // ---------------------------------------
  double dS;      // Step-size for calculation of Jacobian 
  double Smax;    // Maximum step-size
  double Smin;    // Minimum step-size
  double Eps;     // Minimum magnitude of (F1^2+F2^2)^1/2 at root F1(x1,x2) = F2(x1,x2) = 0
  double alpha;   // Ensures sufficient decrease in function value
  int    MaxIter; // Maximum number of iterations

  // ------------------
  // Physical constants
  // ------------------
  double pc_e;         // Magnitude of electron charge (C)
  double pc_c;         // Velocity of light in vacuum (m/s)
  double pc_m_e;       // Electron rest mass (kg)
  double pc_m_p;       // Proton rest mass (kg)
  double pc_epsilon_0; // Vacuum permittivity (F/m)
  double pc_mu_0;      // Vacuum permeability (F/m)

  // ----------------------
  // Mathematical constants
  // ----------------------
  complex<double> mc_i;    // Square root of minus one
  double          mc_spi;  // Square root of pi
  double          mc_sp2;  // Square root of 2 pi
  
  // Constructor
  Utility ();
  // Destructor
  ~Utility ();

  // Evaluate plasma dispersion function
  complex<double> ZPlasma (complex<double> xi);

  // Evaluate right-hand sides of differential equations for RK4/RK5 integration routines
  virtual void CashKarp45Rhs (double x, double* y, double* dydx);
  // Adaptive step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag, 
			   int diag, FILE* file);
  // Fixed step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h);

  // Evaluate right-hand sides of differential equations for RK4/RK5 integration routines
  virtual void CashKarp45Rhs (double x, complex<double>* y, complex<double>* dydx);
  // Advance set of coupled first-order o.d.e.s by single step using adaptive
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Adaptive (int neqns, double& x, complex<double>* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag, 
			   int diag, FILE* file);
  // Advance set of coupled first-order o.d.e.s by single step using fixed
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Fixed (int neqns, double& x, complex<double>* y, complex<double>* err, double h);

  // Evaluate right-hand sides of differential equations for RK4/RK5 integration routines
  virtual void CashKarp45Rhs1 (double x, complex<double>* y, complex<double>* dydx);
  // Advance set of coupled first-order o.d.e.s by single step using adaptive
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Adaptive1 (int neqns, double& x, complex<double>* y, double& h, 
			    double& t_err, double acc, double S, double T, int& rept,
			    int maxrept, double h_min, double h_max, int flag, 
			    int diag, FILE* file);
  // Advance set of coupled first-order o.d.e.s by single step using fixed
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Fixed1 (int neqns, double& x, complex<double>* y, complex<double>* err, double h);

  // Evaluate target function for one-dimensional root finding
  virtual double RootFindF (double x);
  // One-dimensional root finding routine
  double RootFind (double x1, double x2);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double x1, double x2, double F1, double F2, double& x);

  // Evaluate target functions for two-dimensional root finding
  virtual void NewtonFunction (double x1, double x2, double& F1, double& F2);
  // Find root of F1(x1,x2) = F2(x1,x2) = 0 via Newton-Raphson method
  void NewtonRoot (double& x1, double& x2, double& Residual, int verbose);
  // Backtrack along Newton step in order to minimize f = (F1*F1 + F2*F2) /2
  void NewtonBackTrack (double& x1, double& x2, double dx1, double dx2, double& dx, double f, double g1, double g2, double& lambda);
  // Calculate Jacobian matrix for Newton-Raphson root finding
  void NewtonJacobian (double x1, double x2, double& J11, double& J12, double& J21, double& J22);

  // Return maximum of two values
  double Fmax (double f1, double f2);
  // Return minimum of two values
  double Fmin (double f1, double f2);
  
  // Strip comments from a string
  string StripComments (const string& input);
  // Read JSON file
  json ReadJSONFile (const string& filename);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);
  // Open file for appending
  FILE* OpenFilea (char* filename);

  // Check that directory exists, and create it otherwise
  bool CreateDirectory (const char* path);

  // Call operating system
  void CallSystem (char* command);
};
