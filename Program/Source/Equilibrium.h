// Equilibrium.h

// ########################################################################################
// Class to solve inverse aspect-ratio expanded tokamak equilibrium problem.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).
// Radial coordinate, r, normalized to epsa * R_0, where eps_a is inverse-aspect ratio.
// So r=0 is magnetic axis, r=1 is plasma/vacuum interface.

// Flux-surfaces:

// R(r,w) = 1 - epsa r cosw + epsa^2 H1(r) + epsa^2 sum_{n=2,Ns} [Hn(r) cos(n-1)w + Vn(r) sin(n-1)w]
// Z(r,w) =     epsa r sinw                + epsa^2 sum_{n=2,Ns} [Hn(r) sin(n-1)w - Vn(r) cos(n-1)w]
//
// Here, R,phi,Z are cylindrical polar coordinates while r,w,phi are flux coordinates

// Edge shaping: Hna = Hn(1), Vna = Vn(1)

// Equilibrium profiles:
// 
//  f1(r) = (1/nu/qc) [1 - (1 - r^2)^nu]
//
//  p2(r) = pc (1 - r^2)^mu
//
// qc is lowest-order safety-factor on magnetic axis.

// Inputs:
//  Inputs/Namelist.nml - namelist
//  Inputs/Shape.txt    - nlines
//                        H2a V2a
//                        H3a V3a
//                        etc.

// Outputs:
//  Plots/Equilibrium.nc

// Plots:
//  Plots/*.py

// Program uses:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// ########################################################################################

#ifndef EQUILIBRIUM
#define EQUILIBRIUM

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <netcdf>

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;

// Namelist reading function
extern "C" void NameListEquilibrium (double* QC, double* NU, double* PC, double* MU, double* EPSA,
				     double* EPS, int* NS, int* NR, int* NF, int* NW,
				     double* ACC, double* H0, double* HMIN, double* HMAX);
    
// ############
// Class header
// ############
class Equilibrium
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  double qc;     // Lowest-order safety-factor on magnetic axis (read from namelist)
  double nu;     // Toroidal current peaking parameter (read from namelist)
  double pc;     // Normalized plasma pressure on magnetic axis (read from namelist)
  double mu;     // Pressure peaking parameter (read from namelist)
  double epsa;   // Inverse aspect-ratio of plasma (read from namelist)

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from namelist)
  int    Ns;     // Number of shaping harmonics (read from namelist)
  int    Nr;     // Number of radial grid points (read from namelist)
  int    Nf;     // Number of magnetic flux-surfaces (read from namelist)
  int    Nw;     // Number of angular points on magnetic flux-surfaces (read from namelist)

  // ----------------
  // Calculation data
  // ----------------
  double* rr;               // Radial grid points
  double* p2;               // Plasma pressure profile
  double* f1;               // Lowest-order poloidal flux function
  double* f3;               // Higher-order poloidal flux function
  double* g2;               // Toroidal flux function
  double* q0;               // Lowest-order safety-factor
  double* q2;               // Higher-order safety-factor
  double* It;               // Toroidal plasma current
  double* Ip;               // Poloidal plasma current
  double* Jt;               // Radial derivative of toroidal plasma current
  double* Jp;               // Radial derivative of poloidal plasma current
  double* pp;               // Radial plasma pressure gradient
  double* ppp;              // Second radial derivative of plasma pressure 
  double* qq;               // First radial derivative of safety-factor times r
  double* qqq;              // Radial derivative of qq times r
  double* s;                // Magnetic shear:     s  = r q'/q
  double* s2;               // Higher-order shear: s2 = r^2 q''/q
  double* S1;               // First shaping function
  double* S2;               // Second shaping function
  double* P1;               // First profile function: (2-s)/q
  double* P2;               // Second profile function: r dP1/dr
  double* P3;               // Third profile function
  double* P3a;              // Auxillary third profile function
  double* ff;               // f profile
  double* ggr2;             // <|nabla r|^2> profile
  double* RR2;              // <R^2> profile
   
  Array<double,2> HHfunc;   // Horizontal shaping functions
  Array<double,2> VVfunc;   // Vertical shaping functions
  Array<double,2> HPfunc;   // Radial derivatives of horizontal shaping functions
  Array<double,2> VPfunc;   // Radial derivatives of vertical shaping functions

  gsl_spline* Itspline;     // Interpolated It function
  gsl_spline* Ipspline;     // Interpolated Ip function

  gsl_interp_accel* Itacc;  // Accelerator for interpolated It function
  gsl_interp_accel* Ipacc;  // Accelerator for interpolated Ip function

  gsl_spline*  g2spline;    // Interpolated g2 function
  gsl_spline** HHspline;    // Interpolated horizontal shaping functions
  gsl_spline** VVspline;    // Interpolated vertical shaping functions
  gsl_spline** HPspline;    // Interpolated radial derivatives of horizontal shaping functions
  gsl_spline** VPspline;    // Interpolated radial derivatives of vertical shaping functions

  gsl_interp_accel*  g2acc; // Accelerator for interpolated g2 function
  gsl_interp_accel** HHacc; // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc; // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc; // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc; // Accelerator for interpolated radial derivatives of vertical shaping functions

  gsl_spline* fspline;      // Interpolated f function
  gsl_spline* gr2spline;    // Interpolated <|nabla r|^2> function
  gsl_spline* R2spline;     // Interpolated <R^2> function

  gsl_interp_accel* facc;   // Accelerator for interpolated f function
  gsl_interp_accel* gr2acc; // Accelerator for interpolated <|nabla r|^2> function
  gsl_interp_accel* R2acc;  // Accelerator for interpolated <R^2> function

  Array<double,2> RR;       // R coodinates of magnetic flux-surfaces
  Array<double,2> ZZ;       // Z coodinates of magnetic flux-surfaces
   
  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;     // Integration accuracy (read from namelist)
  double h0;      // Initial integration step-length (read from namelist)
  double hmin;    // Minimum integration step-length (read from namelist)
  double hmax;    // Maximum integration step-length (read from namelist)
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Integration error calcualation flag
  
  // ----------------------------
  // Cash-Karp RK4/RK5 parameters
  // ----------------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int count, rhs_chooser;

 public:

  // ..................
  // in Equilibrium.cpp
  // ..................
  
  // Constructor
  Equilibrium ();
  // Destructor
  ~Equilibrium ();

  // Solve problem
  void Solve ();

private:

  // ..................
  // in Equilibrium.cpp
  // ..................
  
  // Return f1(r)
  double Getf1 (double r);
  // Return f1'(r)
  double Getf1p (double r);
  // Return p2(r)
  double Getp2 (double r);
  // Return p2'(r)
  double Getp2p (double r);
  // Return p2''(r)
  double Getp2pp (double r);
  
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, double* y, double* dydx);

  // Adaptive step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag, 
			   int diag, FILE* file);
  // Fixed step-length Cash-Karp RK4/RK5 integration routine
  void CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);
};

#endif //EQUILIBRIUM
