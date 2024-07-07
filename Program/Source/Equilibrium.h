// Equilibrium.h

// ########################################################################################
// Class to solve inverse aspect-ratio expanded tokamak equilibrium problem.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).
// Radial coordinate, r, normalized to epsa * R_0, where eps_a is inverse-aspect ratio.
// So r = 0. is magnetic axis, r = 1. is plasma/vacuum interface.

// Flux-surfaces:

// R(r,w) = 1 - epsa r cosw + epsa^2 H1(r) + epsa^2 sum_{n=2,Ns} [Hn(r) cos(n-1)w + Vn(r) sin(n-1)w]
// Z(r,w) =     epsa r sinw                + epsa^2 sum_{n=2,Ns} [Hn(r) sin(n-1)w - Vn(r) cos(n-1)w]
//
// Here, R, phi, Z are cylindrical polar coordinates while r, w, phi are flux coordinates

// Edge shaping: Hna = Hn(1.), Vna = Vn(1.)

// Equilibrium profiles:
//
//  Lowest order (i.e., cylindrical) safety factor profile is q0(r) = r^2/f1(r)
//  Pressure profile is p(r) = epsa^2 p2(r)
// 
//  f1(r) = (1/nu/qc) [1 - (1 - r^2)^nu]
//
//  p2(r) = pc (1 - r^2)^mu
//
// qc is lowest-order safety-factor on magnetic axis.

// Inputs:
//  Inputs/Namelist.nml - namelist
//  Inputs/Shape.txt    - nlines (i.e. no of subsequent lines to be read)
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
				     double* EPS, int* NS, int* NR, int* NF, int* NW, int* HIGH, 
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
  double epsa;   // Inverse aspect-ratio of plasma (read from namelist)
  double qc;     // Lowest-order safety-factor on magnetic axis (read from namelist)
  double nu;     // Toroidal current peaking parameter (read from namelist)
  double pc;     // Normalized plasma pressure on magnetic axis (read from namelist)
  double mu;     // Pressure peaking parameter (read from namelist)

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps;    // Distance of closest approach to magnetic axis (read from namelist)
  int    Ns;     // Number of shaping harmonics (read from namelist)
  int    Nr;     // Number of radial grid-points for calculation purposes (read from namelist)
  int    Nf;     // Number of radial grid-points for visualization purposes (read from namelist)
  int    Nw;     // Number of angular grid-points for visulalization purposes (read from namelist)
  int    HIGH;   // Switch for higher order calculation of boundary data (read from namelist)
                 //  HIGH = 0 - lowest-order analytic calculation of boundary data
                 //  HIGH = 1 - higher-order analytic calculation of boundary data
                 //  HIGH = 2 - numerical calculation of boundary data

  // ----------------
  // Calculation data
  // ----------------
  double* rr;    // Radial grid-points
  double* p2;    // Plasma pressure profile
  double* f1;    // Lowest-order poloidal flux function
  double* f3;    // Higher-order poloidal flux function
  double* g2;    // Toroidal flux function
  double* q0;    // Lowest-order safety-factor
  double* q2;    // Higher-order safety-factor
  double* It;    // Toroidal plasma current
  double* Ip;    // Poloidal plasma current
  double* Jt;    // Radial derivative of toroidal plasma current
  double* Jp;    // Radial derivative of poloidal plasma current
  double* pp;    // Radial derivative of plasma pressure
  double* ppp;   // Second radial derivative of plasma pressure 
  double* qq;    // First radial derivative of safety-factor times r
  double* qqq;   // Radial derivative of qq times r
  double* s;     // Magnetic shear:              s  = r q2'/q2
  double* s2;    // Second-order magnetic shear: s2 = r^2 q2''/q2
  double* S1;    // First shaping function
  double* S2;    // Second shaping function
  double* P1;    // First profile function:  (2-s)/q2
  double* P2;    // Second profile function: r dP1/dr
  double* P3;    // Third profile function
  double* P3a;   // Auxillary third profile function
  double* ff;    // f profile
  double* ggr2;  // <|nabla r|^2> profile
  double* RR2;   // <R^2> profile
   
  Array<double,2>    HHfunc;    // Horizontal shaping functions
  Array<double,2>    VVfunc;    // Vertical shaping functions
  Array<double,2>    HPfunc;    // Radial derivatives of horizontal shaping functions
  Array<double,2>    VPfunc;    // Radial derivatives of vertical shaping functions
  double*            Lfunc;     // Relabelling function

  gsl_spline*        Itspline;  // Interpolated It function
  gsl_spline*        Ipspline;  // Interpolated Ip function

  gsl_interp_accel*  Itacc;     // Accelerator for interpolated It function
  gsl_interp_accel*  Ipacc;     // Accelerator for interpolated Ip function

  gsl_spline*        g2spline;  // Interpolated g2 function
  gsl_spline**       HHspline;  // Interpolated horizontal shaping functions
  gsl_spline**       VVspline;  // Interpolated vertical shaping functions
  gsl_spline**       HPspline;  // Interpolated radial derivatives of horizontal shaping functions
  gsl_spline**       VPspline;  // Interpolated radial derivatives of vertical shaping functions
  gsl_spline*        Lspline;   // Interpolated relabelling function
  gsl_spline*        wspline;   // Interpolated omega function

  gsl_interp_accel*  g2acc;     // Accelerator for interpolated g2 function
  gsl_interp_accel** HHacc;     // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc;     // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc;     // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc;     // Accelerator for interpolated radial derivatives of vertical shaping functions
  gsl_interp_accel*  Lacc;      // Accelerator for interpolated relabelling function
  gsl_interp_accel*  wacc;      // Accelerator for interpolated omega function

  gsl_spline*        fspline;   // Interpolated f function
  gsl_spline*        gr2spline; // Interpolated <|nabla r|^2> function
  gsl_spline*        R2spline;  // Interpolated <R^2> function

  gsl_interp_accel*  facc;      // Accelerator for interpolated f function
  gsl_interp_accel*  gr2acc;    // Accelerator for interpolated <|nabla r|^2> function
  gsl_interp_accel*  R2acc;     // Accelerator for interpolated <R^2> function

  Array<double,2>    RR;        // R coodinates of magnetic flux-surfaces for visualization purposes
  Array<double,2>    ZZ;        // Z coodinates of magnetic flux-surfaces for visualization purposes
  Array<double,2>    rvals;     // r values on magnetic flux-surfaces for visualization purposes
  Array<double,2>    thvals;    // theta values on magnetic flux-surfaces for visualization purposes

  double*            Rbound;    // R values on plasma boundary
  double*            Zbound;    // Z values on plasma boundary
  double*            tbound;    // theta values on plasma boundary
  double*            wbound;    // omega values on plasma boundary
  double*            tbound0;   // Preliminary theta values on plasma boundary
  double*            wbound0;   // Preliminary omega values on plasma boundary
  double*            R2b;       // R^2 values on plasma boundary
  double*            grr2b;     // |nabla r|^2 values on plasma boundary
   
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

  // Return relabelling parameter
  double GetL (double r);
  // Return R
  double GetR (double r, double w);
  // Return dRdr
  double GetdRdr (double r, double w);
  // Return dRdw
  double GetdRdw (double r, double w);
  // Return Z
  double GetZ (double r, double w);
  // Return dZdr
  double GetdZdr (double r, double w);
  // Return dZdw
  double GetdZdw (double r, double w);
  // Return w-theta transformation function
  double Gettfun (double r, double w);
  // Return R2
  double GetR2 (double r, double t);
  // Return |nabla r|^2
  double Getgrr2 (double r, double t);
  
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
