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
//  nclohmann JSON library (https://github.com/nlohmann/json)

// Source: https://github.com/rfitzp/TJ

// ############################################################################################

#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

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

using namespace std;
using           json = nlohmann::json;

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

  // -----------------------
  // Root finding parameters
  // -----------------------
  int    Nint;    // Number of search intervals
  double Eta;     // Min. magnitude of f at root f(x) = 0
  int    Maxiter; // Maximum number of iterations
  
  // ..............
  // in Utility.cpp
  // ..............
  
  // Constructor
  Utility ();
  // Destructor
  ~Utility ();

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

  // Target function for one-dimensional root finding
  virtual double RootFindF (double x);
  // One-dimensional root finding routine
  double RootFind (double x1, double x2);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double x1, double x2, double F1, double F2, double& x);

  // Strip comments from a string
  string StripComments (const string& input);
  // Read JSON file
  json ReadJSONFile (const string& filename);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);

  // Check that directory exists, and create it otherwise
  bool CreateDirectory (const char* path);

  // Call operating system
  void CallSystem (char* command);
};
