// TJ.h

// #################################################################################
// Class to calculate tearing stability matrix and tearing eigenfunctions in an
// inverse aspect-ratio expanded tokamak equilibrium.

// See Documentation/TJ.tex

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field).
// Radial coordinate, r, normalized to epsa * R_0, where eps_a is inverse-aspect ratio.
// So r=0 corresponds to magnetic axis, and r=1 to plasma/vacuum interface.

// Program assumes monotonic safety-factor profile.

// Inputs:
//  Input/Namelist.nml - namelist

// Outputs:
//  Plots/TJ.nc

// Plots:
//  Plots/*.py

// Program uses:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)
//  Armadillo library      (https://arma.sourceforge.net)

// #################################################################################

#ifndef TJXX
#define TJXX

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <netcdf>
#include <armadillo>

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace arma;

// Namelist reading function
extern "C" void NameListTJ (int* NTOR, int* MMIN, int* MMAX,
			    double* EPS, double* DEL, int* NFIX, int* NDIAG, double* NULC, int* ITERMAX, int* FREE, 
			    double* ACC, double* H0, double* HMIN, double* HMAX, double* EPSF);

// ############
// Class header
// ############
class TJ
{
 private:

  // ..................
  // Control parameters
  // ..................
  int    NTOR;    // Toroidal mode number of RMP (read from namelist)
  int    MMIN;    // Minimum poloidal mode number included in calculation (read from namelist)
  int    MMAX;    // Maximum poloidal mode number included in calculation (read from namelist)

  double EPS;     // Solutions launched from magnetic axis at r = EPS (read from namelist)
  double DEL;     // Distance of closest approach to rational surface is DEL (read from namelist)
  int    NFIX;    // Number of fixups (read from namelist)
  int    NDIAG;   // Number of radial grid-points for diagnostics (read from namelist)
  double NULC;    // Use zero pressure jump conditions when |nu_L| < NULC (read from namelist)
  int    ITERMAX; // Maximum number of iterations used to determine quantities at rational surface (read from namelist)
  int    FREE;    // Flag for free/fixed boundary calculation 

  double EPSF;    // Step-length for finite difference determination of derivative

  // ..................................................
  // Equilibrium data (read from Inputs/Equilibrium.nc)
  // ..................................................
  double             epsa;      // Inverse aspect-ratio
  int                Ns;        // Number of shaping harmonics
  int                Nr;        // Number of radial grid points

  double*            rr;        // Radial grid points
  double*            pp;        // First radial derivative of plasma pressure
  double*            ppp;       // Second radial derivative of plasma pressure 
  double*            q;         // Safety-factor
  double*            s;         // Magnetic shear
  double*            s2;        // Higher-order shear: s2 = r^2 q''/q
  double*            S1;        // First shaping function
  double*            P1;        // First profile function: (2-s)/q
  double*            P2;        // Second profile function: r dP1/dr
  double*            P3;        // Third profile function

  Array<double,2>    HHfunc;    // Horizontal shaping functions
  Array<double,2>    VVfunc;    // Vertical shaping functions
  Array<double,2>    HPfunc;    // Radial derivatives of horizontal shaping functions
  Array<double,2>    VPfunc;    // Radial derivatives of vertical shaping functions

  gsl_spline*        ppspline;  // Interpolated pp function
  gsl_spline*        pppspline; // Interpolated ppp function
  gsl_spline*        qspline;   // Interpolated q function
  gsl_spline*        sspline;   // Interpolated s function
  gsl_spline*        s2spline;  // Interpolated s2 function
  gsl_spline*        S1spline;  // Interpolated S1 function
  gsl_spline*        P1spline;  // Interpolated P1 function
  gsl_spline*        P2spline;  // Interpolated P2 function
  gsl_spline*        P3spline;  // Interpolated P3 function
  
  gsl_spline**       HHspline;  // Interpolated horizontal shaping functions
  gsl_spline**       VVspline;  // Interpolated vertical shaping functions
  gsl_spline**       HPspline;  // Interpolated radial derivatives of horizontal shaping functions
  gsl_spline**       VPspline;  // Interpolated radial derivatives of vertical shaping functions

  gsl_interp_accel*  ppacc;     // Accelerator for interpolated pp function
  gsl_interp_accel*  pppacc;    // Accelerator for interpolated ppp function
  gsl_interp_accel*  qacc;      // Accelerator for interpolated q function
  gsl_interp_accel*  sacc;      // Accelerator for interpolated s function
  gsl_interp_accel*  s2acc;     // Accelerator for interpolated s2 function
  gsl_interp_accel*  S1acc;     // Accelerator for interpolated S1 function
  gsl_interp_accel*  P1acc;     // Accelerator for interpolated P1 function
  gsl_interp_accel*  P2acc;     // Accelerator for interpolated P2 function
  gsl_interp_accel*  P3acc;     // Accelerator for interpolated P3 function
   
  gsl_interp_accel** HHacc;     // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc;     // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc;     // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc;     // Accelerator for interpolated radial derivatives of vertical shaping functions

  int             Nf;           // Number of magnetic flux-surfaces
  int             Nw;           // Number of angular points on magnetic flux-surfaces
  Array<double,2> RR;           // R coodinates of magnetic flux-surfaces
  Array<double,2> ZZ;           // Z coodinates of magnetic flux-surfaces
  Array<double,2> rvals;        // r values on magnetic flux-surfaces
  Array<double,2> thvals;       // theta values on magnetic flux-surfaces

  // ......................
  // Calculation parameters
  // ......................
  double   ntor;  // Toroidal mode number
  int      J;     // Number of poloidal harmonics included in calculation
  int      K;     // Number of solution vectors: K = J + nres
  int*     MPOL;  // Poloidal mode numbers of included poloidal harmonics
  double*  mpol;  // Poloidal mode numbers of included poloidal harmonics

  // .....................
  // Rational surface data
  // ....................
  int     nres;    // Number of rational magnetic flux-surfaces
  int*    mres;    // Poloidal mode numbers at rational surfaces
  double* qres;    // Safety-factors at rational surfaces
  double* rres;    // Minor radii of rational surfaces
  double* qerr;    // Residual in determination of rational surface
  double* sres;    // Magnetic shears at rational surfaces
  double* DIres;   // DI values at rational surfaces
  double* nuLres;  // Mercier indices of large solution at rational surfaces
  double* nuSres;  // Mercies indices of small solution at rational surfaces
  int*    Jres;    // Index of resonant poloidal harmonic at rational surfaces

  // --------------------
  // Vacuum solution data
  // --------------------
  double                   sa;    // Edge magnetic shear
  double                   G1;    // Vacuum solution parameter
  double                   G2;    // Vacuum solution parameter
  Array<complex<double>,2> Pvac;  // Vacuum solution matrix
  Array<complex<double>,2> Qvac;  // Vacuum solution matrix
  Array<complex<double>,2> Rvac;  // Vacuum solution matrix
  Array<complex<double>,2> Svac;  // Vacuum solution matrix
  Array<complex<double>,2> Pdag;  // Hermitian conjugate of Pvac
  Array<complex<double>,2> Rdag;  // Hermitian conjugate of Rvac
  Array<complex<double>,2> Avac;  // Vacuum residual matrix
  Array<complex<double>,2> Bvac;  // Vacuum residual matrix
  Array<complex<double>,2> Cvac;  // Vacuum residual matrix
  Array<complex<double>,2> Hmat;  // Vacuum response matrix
  Array<complex<double>,2> Hdag;  // Hermitian conjugate of Hmat
  Array<complex<double>,2> Hsym;  // Symmeterized Hmat

  // -----------------
  // ODE Solution data
  // -----------------
  double*                  Rgrid;  // Radial grid points for diagnostics
  Array<double,2>          Ttest;  // Torque test for solution vectors versus radius
  Array<double,2>          Pnorm;  // Norms of Psi components of solution vectors versus radius
  Array<double,2>          Znorm;  // Norms of Z components of solution vectors versus radius
  double*                  hode;   // Step-length versus radius
  double*                  eode;   // Truncation error versus radius
  Array<complex<double>,3> YYY;    // Solution vectors versus radius
  Array<complex<double>,2> Pi;     // Reconnected fluxes at rational surfaces associated with solution vectors

  // -------------------------------------
  // Tearing mode dispersion relation data
  // -------------------------------------
  Array<complex<double>,2> Psia;  // Values of Psi at plasma boundary due to continuous solutions launched from magnetic axis
  Array<complex<double>,2> Za;    // Values of Z at plasma boundary due to continuous solutions launched from magnetic axis
  Array<complex<double>,2> Pia;   // Values of reconnected fluxes at rational surfaces due to continuous solutions launched from magnetic axis
  Array<complex<double>,2> Psis;  // Values of Psi at plasma boundary due to small solutions launched from rational surfaces
  Array<complex<double>,2> Zs;    // Values of Z at plasma boundary due to small solutions launched from rational surfaces
  Array<complex<double>,2> Pis;   // Values of reconnected fluxes at rational surfaces due to small solutions launched from rational surfaces
  Array<complex<double>,2> Xmat;  // X-matrix
  Array<complex<double>,2> Ymat;  // Y-matrix
  Array<complex<double>,2> Omat;  // Omega-matrix
  Array<complex<double>,2> Fmat;  // Inductance matrix
  Array<complex<double>,3> Psif;  // Psi components of fully reconnected tearing eigenfunctions
  Array<complex<double>,3> Zf;    // Z componnents of fully reconnected tearing eigenfunctions
  Array<complex<double>,2> Emat;  // Tearing stability matrix
  Array<complex<double>,3> Psiu;  // Psi components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> Zu;    // Z componnents of unreconnected tearing eigenfunctions
  Array<double,2>          Tf;    // Torques associated with fully reconnected eigenfunctions
  Array<double,2>          Tu;    // Torques associated with unreconnected eigenfunctions
  Array<double,3>          Tfull; // Torques associated with pairs of fully reconnected eigenfunctions
  Array<double,3>          Tunrc; // Torques associated with pairs of unreconnected eigenfunctions
    
  // .......................
  // Root finding parameters
  // .......................
  double Eta;      // Minimum magnitude of f at root f(x) = 0
  int    Maxiter;  // Maximum number of iterations

  // ...............................
  // Adaptive integration parameters
  // ...............................
  double acc;     // Integration accuracy (read from namelist)
  double h0;      // Initial step-length (read from namelist)
  double hmin;    // Minimum step-length (read from namelist)
  double hmax;    // Maximum step-length (read from namelist)
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
  int    count;
  double qval;
  
 public:

  // .........
  // In TJ.cpp
  // .........

  // Constructor
  TJ ();
  // Destructor
  ~TJ ();

  // Solve problem
  void Solve ();

 private:

  // .........
  // In TJ.cpp
  // .........
  
  // Set mode numbers
  void SetModeNumbers ();
  // Write data to Plots/TJ.nc
  void WriteNetCDF ();
  // Deallocate memory
  void CleanUp ();
  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);

  // ..................
  // In Equilibrium.cpp
  // ..................

  // Read equilibrium data from Inputs/Equilibrium.nc
  void ReadEquilibrium ();
    
  // .............
  // In Vacuum.cpp
  // .............

  // Calculate vacuum matrices
  void GetVacuum ();
  
  // ...............
  // In Rational.cpp
  // ...............

  // Target function for finding rational surfaces
  double Feval (double r);
  // Find rational surfaces
  void FindRational ();

  // .............
  // In Matrix.cpp
  // .............

  // Get values of coupling matrices
  void GetMatrices (double r,
		    Array<complex<double>,2> LLmmp, Array<complex<double>,2> MMmmp,
		    Array<complex<double>,2> NNmmp, Array<complex<double>,2> PPmmp);
  // Get values of coupling matrix elements
  void GetMatrices (double r, int m, int mp,
		    complex<double>& Lmmp, complex<double>& Mmmp,
		    complex<double>& Nmmp, complex<double>& Pmmp);

  // ...............
  // In Resonant.cpp
  // ...............
  
  // Calculate coefficients L1, P1, and T1 at given radius
  void GetL1P1T1 (double rm, double mm, int km, double& sm, double& L1, double& P1, double& T1);
  // Calculate coefficients L1k, M1k, N1k, and P1k at given radius
  void GetL1kM1kN1kP1k (double rm, double mm, int km, complex<double>* L1k, complex<double>* M1k,
			complex<double>* N1k, complex<double>* P1k);
  // Jump solution vectors across jresth rational surface
  void JumpRational (int jres, double& r, Array<complex<double>,2> YY);
  // Calculate angular momentum flux associated with ith solution vector
  double GetTorque (double r, Array<complex<double>,2> Psi, Array<complex<double>,2> Z, int i);
  
  // ...............
  // In ODESolve.cpp
  // ...............

  // Function to solve outer region odes
  void ODESolve ();
  // Launch solution vectors from magnetic axis
  void LaunchAxis (double r, Array<complex<double>,2> YY);
  // Launch solution vector from jresth rational surface
  void LaunchRational (int jres, double r, Array<complex<double>,2> YY);
  // Integrate multiple solution vectors from given value of r to r = rx while performing nf fixups
  void SegmentFixup (double& r, double rx, int nf, Array<complex<double>,2> YY);
  // Integrate multiple solution vectors from given value of r to to r = rx
  void Segment (double& r, double rx, Array<complex<double>,2> YY);
  // Perform fixup of multiple solution vectors lauched from magnetic axis
  void Fixup (double r, Array<complex<double>,2> YY);
  // Evaluate right-hand sides of outer region odes
  void Rhs (double r, complex<double>* Y, complex<double>* dYdr);
  // Pack YY solution vector
  void PackYY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, Array<complex<double>,2> YY);
  // Unpack YY solution vector
  void UnpackYY (Array<complex<double>,2> YY, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ);
  // Pack YY solution vectors 
  void PackYY (complex<double>* Y, Array<complex<double>,2> YY);
  // Unpack YY solution vector
  void UnpackYY (Array<complex<double>,2> YY, complex<double>* Y);
  // Pack Y solution vectors
  void PackY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, complex<double>* Y);
  // Unpack Y solution vectors
  void UnpackY (complex<double>* Y, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ);
  // Function to perform torque test on solution vectors
  double TorqueTest (double r, Array<complex<double>,2> YY);
  // Calculate angular momentum flux associated with ith solution vector
  double GetTorque (double r, Array<complex<double>,2> YY, int i);
  // Calculate norms of ith solution vector
  void GetNorms (Array<complex<double>,2> YY, int i, double& Pnorm, double &Znorm);
 
  // .................
  // In Dispersion.cpp
  // .................
  
  // Find tearing mode dispersion relation and construct tearing eigenfunctions
  void FindDispersion ();
  // Calculate angular momentum flux associated with pairs of fully reconnected solution vectors
  void GetTorqueFull ();
  // Calculate angular momentum flux associated with pairs of unreconnected solution vectors
  void GetTorqueUnrc ();
  
  // ..................
  // In Interpolate.cpp
  // ..................

  // Return value of pp
  double Getpp (double r);
  // Return value of ppp
  double Getppp (double r);
  // Return value of q
  double Getq (double r);
  // Return value of s
  double Gets (double r);
  // Return value of s2
  double Gets2 (double r);
  // Return value of S1
  double GetS1 (double r);
  // Return value of P1
  double GetP1 (double r);
  // Return value of P2
  double GetP2 (double r);
  // Return value of P3
  double GetP3 (double r);
  // Return value of Hn
  double GetHn (int n, double r);
  // Return value of Hnp
  double GetHnp (int n, double r);
  // Return value of Vn
  double GetVn (int n, double r);
  // Return value of Vnp
  double GetVnp (int n, double r);
  // Return value of DI
  double GetDI (double r);

  // ...............
  // In ZeroFind.cpp
  // ...............

  //  Routine to find approximate root of F(x) = 0 using Ridder's method. 
  double RootFind ();
  // Ridder's method for finding root of F(x) = 0.
  void Ridder (double x1, double x2, double F1, double F2, double& x);
  
  // ................
  // In Integrate.cpp
  // ................

  // Advance set of coupled first-order o.d.e.s by single step using adaptive
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Adaptive (int neqns, double& x, complex<double>* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag, 
			   int diag, FILE* file);
  // Advance set of coupled first-order o.d.e.s by single step using fixed
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Fixed (int neqns, double& x, complex<double>* y, complex<double>* err, double h);
  
  // ................
  // In Armadillo.cpp
  // ...............

  // Solve linear system of equations A . X = B, for B, where all quantities are complex rectangular matrices
  void SolveLinearSystem (Array<complex<double>,2> A, Array<complex<double>,2> X, Array<complex<double>,2> B);
  // Invert square matrix
  void InvertMatrix (Array<complex<double>,2> A, Array<complex<double>,2> invA);
};

#endif //TJXX
