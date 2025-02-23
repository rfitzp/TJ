// LightEquilibrium.h

// ########################################################################################

// Lightweight class to solve inverse aspect-ratio expanded tokamak equilibrium problem.
// Main class method determines nu value that corresponds to given edge safety-factor.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).
// Radial coordinate, r, normalized to epsa * R_0, where epsa is inverse-aspect ratio.
// So r = 0 is magnetic axis and r = 1 is plasma/vacuum interface.

// Equilibrium magnetic flux-surfaces:

// R(r,w) = 1 - epsa r cosw + epsa^2 H1(r) + epsa^2 sum_{n=2,Ns} [Hn(r) cos(n-1)w + Vn(r) sin(n-1)w]
// Z(r,w) =     epsa r sinw                + epsa^2 sum_{n=2,Ns} [Hn(r) sin(n-1)w - Vn(r) cos(n-1)w]
//
// Here, R, phi, Z are cylindrical polar coordinates, epsa is the inverse aspect-ratio, r is a flux-surface label,
// and w is a poloidal angle. The class also uses the r, theta, phi (PEST) straight field-line coordinate system
// whose Jacobian is r R^2. 

// Edge shaping: Hna = Hn(1), Vna = Vn(1), etc.

// Equilibrium profiles:

// Lowest order (i.e., cylindrical) safety factor profile is q0(r) = r^2 /f1(r)
// Pressure profile is P(r) = epsa^2 p2(r)
// 
//  f1(r) = (1 /nu/qc) [1 - (1 - r^2)^nu] 
//
//  p2(r) = pc (1 - r^2)^mu
//
// qc is lowest-order safety-factor on magnetic axis.
// nu * qc is lowest-order safety-factor at plasma/vacuum interface.

// Inputs:
//  Inputs/Equilibrium.json - JSON file

// Class uses following external libraries:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics
//  University of Texas at Austin
//  rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/TJPaper/TJ.pdf

// ########################################################################################

#pragma once

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include "Utility.h"

using namespace blitz;
    
// ############
// Class header
// ############
class LightEquilibrium : private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  double epsa;        // Inverse aspect-ratio of plasma (passed from class Equilibrium)
  double qc;          // Lowest-order safety-factor on magnetic axis (passed from class Equilibrium)
  double pc;          // Normalized plasma pressure on magnetic axis (passed from class Equilibrium)
  double mu;          // Pressure peaking parameter (read from JSON file)
  vector<double> Hna; // H2(1), H3(1), etc (passed from class Equilibrium)
  vector<double> Vna; // V2(1), V3(1), etc (passed from class Equilibrium)

  double qa;          // Target edge safety-factor value (passed from class Equilibrium)
  double nu;          // Toroidal current peaking parameter

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from JSON file)
  int    Ns;     // Number of shaping harmonics (read from JSON file)
  int    Nr;     // Number of radial grid-points for calculation purposes (read from JSON file)
  
  // ----------------
  // Calculation data
  // ----------------
  double* rr;    // Radial grid-points
  double* f1;    // Lowest-order poloidal flux function
  double* f3;    // Higher-order poloidal flux function
  double* g2;    // Toroidal flux function
  double* q0;    // Lowest-order safety-factor
  double* q2;    // Higher-order safety-factor
   
  Array<double,2>    HHfunc;    // Horizontal shaping functions
  Array<double,2>    VVfunc;    // Vertical shaping functions
  Array<double,2>    HPfunc;    // Radial derivatives of horizontal shaping functions
  Array<double,2>    VPfunc;    // Radial derivatives of vertical shaping functions

  gsl_spline*        g2spline;  // Interpolated g2(r) function
  gsl_spline**       HHspline;  // Interpolated horizontal shaping functions versus r
  gsl_spline**       VVspline;  // Interpolated vertical shaping functions versus r
  gsl_spline**       HPspline;  // Interpolated radial derivatives of horizontal shaping functions versus r
  gsl_spline**       VPspline;  // Interpolated radial derivatives of vertical shaping functions versus r

  gsl_interp_accel*  g2acc;     // Accelerator for interpolated g2 function
  gsl_interp_accel** HHacc;     // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc;     // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc;     // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc;     // Accelerator for interpolated radial derivatives of vertical shaping functions

  // ----
  // Misc
  // ----
  int rhs_chooser, fun_chooser;

 public:

  // .......................
  // in LightEquilibrium.cpp
  // .......................
  
  // Constructor
  LightEquilibrium (double qc, double epsa, double pc, vector<double>& Hna, vector<double>& Vn);
  // Destructor
  ~LightEquilibrium ();

  // Calculate nu value that gives required edge safety-factor
  void GetNu (double qa, double& nu_);

private:

  // .......................
  // in LightEquilibrium.cpp
  // .......................

  // Function to return central and edge q-values
  void GetSafety (double _nu, double& qcentral, double& qedge, double& sa, double &sat);
  
  // Return f1(r)
  double Getf1 (double r);
  // Return f1'(r)
  double Getf1p (double r);
  // Return p2'(r)
  double Getp2p (double r);

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;

  // Target function for 1-dimensional root finding
  double RootFindF (double) override;
};
