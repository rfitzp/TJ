// Vertical.h

// ###################################################################################

// Class to calculate vertical stability eigenfunctions in an inverse aspect-ratio
// expanded tokamak equilibrium.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field).
// Radial coordinate, r, normalized to epsa * R_0, where epsa is inverse-aspect ratio.
// So r = 0 corresponds to magnetic axis, and r = 1 to plasma/vacuum interface.

// See Equilibrium.h for description of equilibrium.
// Program assumes monotonic safety-factor profile.

// Inputs:
//  Inputs/Equilibrium.json  - JSON file
//  Inputs/TJ.json           - JSON file
//  Inputs/Vertical.json     - JSON file

// Outputs:
//  Outputs/Vertical/Vertical.nc

// Plotting scripts:
//  Plots/Vertical/*.py

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)
//  Armadillo library      (https://arma.sourceforge.net)

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics,
//  University of Texas at Austin,
//  rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ/

// Documentation: ../Documentation/Vertical.pdf

// ###################################################################################

#pragma once

#define ARMA_WARN_LEVEL 0

#include <time.h>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <netcdf>
#include <armadillo>

#include "Utility.h"
#include "Equilibrium.h"

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace arma;

// ############
// Class header
// ############
class Vertical : private Utility
{
 private:

  // -----------------
  // Calculation flags
  // -----------------
  int SRC;     // Flag for reading profile data from file (read from Equilibrium JSON file)
  int EQLB;    // Flag for equilibrium calculation only (read from Vertical JSON file)
  int VIZ;     // Flag for generating eigenfunction visualization data (read from TJ JSON file)

  // ----------------------
  // Calculation parameters
  // ----------------------
  int            MMIN;    // Minimum poloidal mode number included in calculation (read from Vertical JSON file)
  int            MMAX;    // Maximum poloidal mode number included in calculation (read from Vertical JSON file)

  double         EPS;     // Solutions launched from magnetic axis at r = EPS (read from Vertical JSON file)
  int            NFIX;    // Number of fixups (read from Vertical JSON file)
  int            NDIAG;   // Number of radial grid-points for diagnostics (read from Vertical JSON file)

  double         EPSF;    // Step-length for finite difference determination of derivatives

  // ------------------
  // Machine parameters
  // ------------------
  double B0;      // On-axis toroidal magnetic field-strength (T) (read from Equilibrium JSON file)
  double R0;      // On-axis plasma major radius (m) (read from  Equilibrium JSON file)
  double n0;      // On-axis electron number density (m^-3) (read from  Equilibrium JSON file)
  double alpha;   // Assumed electron number density profile: n0 (1 - r^2)^alpha (read from Equilibrium JSON file)
  double Zeff;    // Effective ion charge number (read from Equilibrium JSON file)
  double Mion;    // Ion mass number (read from Equilibrium JSON file)
  double Chip;    // Perpendicular momentum/energy diffusivity (m^2/s) (read from Equilibrium JSON file)
  double Teped;   // Electron temperature at edge of plasma (eV) (read from Equilibrium JSON file)
  double neped;   // Electron number density at edge of plasma (m^-3) (read from Equilibrium JSON file)
  double apol;    // Plasma minor radius (m)

  // ----------------
  // Mode number data
  // ----------------
  int      J;    // Number of poloidal harmonics included in calculation
  int*     MPOL; // Poloidal mode numbers of included poloidal harmonics
  double*  mpol; // Poloidal mode numbers of included poloidal harmonics

  // ---------------------------------------------------------------
  // Equilibrium data (read from Outputs/Equilibrium/Equilibrium.nc)
  // ---------------------------------------------------------------
  double             epsa;      // Inverse aspect-ratio of plasma
  int                Ns;        // Number of shaping harmonics
  int                Nr;        // Number of radial grid-points

  double*            rr;        // Radial grid-points
  double*            Psi;       // Psi values at grid-points
  double*            PsiN;      // PsiN values at grid-points
  double*            f;         // Radial derivative of poloidal flux 
  double*            g2;        // Second-order toroidal flux
  double*            p2;        // Second-order plasma pressure
  double*            pp;        // First radial derivative of second-order plasma pressure
  double*            ppp;       // Second radial derivative of second-order plasma pressure 
  double*            q;         // Safety-factor
  double*            s;         // Magnetic shear
  double*            s2;        // Higher-order shear: s2 = r^2 q''/q
  double*            s0;        // Lower order magnetic shear
  double*            S1;        // First shaping function
  double*            S2;        // Second shaping function
  double*            S3;        // Third shaping function
  double*            S4;        // Fourth shaping function
  double*            S5;        // Fifth shaping function
  double*            Sig;       // Sigma shaping function
  double*            P1;        // First profile function: (2-s)/q
  double*            P2;        // Second profile function: r dP1/dr
  double*            P3;        // Third profile function
  double*            ne;        // Electron number density
  double*            Te;        // Electron temperature
  double*            nep;       // Radial derivative of electron number density
  double*            Tep;       // Radial derivative of electron temperature

  Array<double,2>    HHfunc;    // Horizontal shaping functions
  Array<double,2>    VVfunc;    // Vertical shaping functions
  Array<double,2>    HPfunc;    // Radial derivatives of horizontal shaping functions
  Array<double,2>    VPfunc;    // Radial derivatives of vertical shaping functions

  gsl_spline*        Pspline;   // Interpolated PsiN function
  gsl_spline*        fspline;   // Interpolated f function
  gsl_spline*        g2spline;  // Interpolated g2 function
  gsl_spline*        p2spline;  // Interpolated p2 function
  gsl_spline*        ppspline;  // Interpolated pp function
  gsl_spline*        pppspline; // Interpolated ppp function
  gsl_spline*        qspline;   // Interpolated q function
  gsl_spline*        sspline;   // Interpolated s function
  gsl_spline*        s2spline;  // Interpolated s2 function
  gsl_spline*        s0spline;  // Interpolated s0 function
  gsl_spline*        S1spline;  // Interpolated S1 function
  gsl_spline*        S2spline;  // Interpolated S2 function
  gsl_spline*        S3spline;  // Interpolated S3 function
  gsl_spline*        S4spline;  // Interpolated S4 function
  gsl_spline*        S5spline;  // Interpolated S5 function
  gsl_spline*        Sigspline; // Interpolated Sigma function
  gsl_spline*        P1spline;  // Interpolated P1 function
  gsl_spline*        P2spline;  // Interpolated P2 function
  gsl_spline*        P3spline;  // Interpolated P3 function
  gsl_spline*        nespline;  // Interpolated ne function
  gsl_spline*        Tespline;  // Interpolated Te function
  gsl_spline*        nepspline; // Interpolated nep function
  gsl_spline*        Tepspline; // Interpolated Tep function

  gsl_interp_accel*  Pacc;      // Accelerator for interpolated P function
  gsl_interp_accel*  facc;      // Accelerator for interpolated f function
  gsl_interp_accel*  g2acc;     // Accelerator for interpolated g2 function
  gsl_interp_accel*  p2acc;     // Accelerator for interpolated p2 function
  gsl_interp_accel*  ppacc;     // Accelerator for interpolated pp function
  gsl_interp_accel*  pppacc;    // Accelerator for interpolated ppp function
  gsl_interp_accel*  qacc;      // Accelerator for interpolated q function
  gsl_interp_accel*  sacc;      // Accelerator for interpolated s function
  gsl_interp_accel*  s2acc;     // Accelerator for interpolated s2 function
  gsl_interp_accel*  s0acc;     // Accelerator for interpolated s0 function
  gsl_interp_accel*  S1acc;     // Accelerator for interpolated S1 function
  gsl_interp_accel*  S2acc;     // Accelerator for interpolated S2 function
  gsl_interp_accel*  S3acc;     // Accelerator for interpolated S3 function
  gsl_interp_accel*  S4acc;     // Accelerator for interpolated S4 function
  gsl_interp_accel*  S5acc;     // Accelerator for interpolated S5 function
  gsl_interp_accel*  Sigacc;    // Accelerator for interpolated Sig function
  gsl_interp_accel*  P1acc;     // Accelerator for interpolated P1 function
  gsl_interp_accel*  P2acc;     // Accelerator for interpolated P2 function
  gsl_interp_accel*  P3acc;     // Accelerator for interpolated P3 function
  gsl_interp_accel*  neacc;     // Accelerator for interpolated ne function
  gsl_interp_accel*  Teacc;     // Accelerator for interpolated Te function
  gsl_interp_accel*  nepacc;    // Accelerator for interpolated nep function
  gsl_interp_accel*  Tepacc;    // Accelerator for interpolated Tep function
  
  gsl_spline**       HHspline;  // Interpolated horizontal shaping functions
  gsl_spline**       VVspline;  // Interpolated vertical shaping functions
  gsl_spline**       HPspline;  // Interpolated radial derivatives of horizontal shaping functions
  gsl_spline**       VPspline;  // Interpolated radial derivatives of vertical shaping functions
   
  gsl_interp_accel** HHacc;     // Accelerator for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc;     // Accelerator for interpolated vertical shaping functions
  gsl_interp_accel** HPacc;     // Accelerator for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc;     // Accelerator for interpolated radial derivatives of vertical shaping functions

  double*            cmu;       // cosh(mu) on plasma boundary
  double*            eeta;      // eta on plasma boundary
  double*            ceta;      // cos(eta) on plasma boundary
  double*            seta;      // sin(eta) on plasma boundary
  double*            R2grgz;    // R^2 nabla r . nabla z   on plasma boundary
  double*            R2grge;    // R^2 nabla r . nabla eta on plasma boundary

  gsl_spline*        Rrzspline; // Interpolated R2grgz function on plasma boundary
  gsl_spline*        Rrespline; // Interpolated R2grge function on plasma boundary
  gsl_spline*        Rbspline;  // Interpolated R function on plasma boundary
  gsl_spline*        Zbspline;  // Interpolated Z function on plasma boundary

  gsl_interp_accel*  Rrzacc;    // Accelerator for interpolated R2grgz function on plasma boundary
  gsl_interp_accel*  Rreacc;    // Accelerator for interpolated R2grge function on plasma boundary
  gsl_interp_accel*  Rbacc;     // Accelerator for interpolated R function on plasma boundary
  gsl_interp_accel*  Zbacc;     // Accelerator for interpolated Z function on plasma boundary

  // --------------------------
  // Plasma boundary parameters
  // --------------------------
  double* tbound; // theta values on plasma boundary
  double* Rbound; // R values on plasma boundary
  double* Zbound; // Z values on plasma boundary
  double* dRdthe; // dR/dtheta values on plasma boundary
  double* dZdthe; // dZ/dtheta values on plasma boundary

  // -------------------------------------
  // Visualization of ideal eigenfunctions
  // -------------------------------------
  int             Nf;     // Number of radial grid-points on visualization grid (from Equilibrium.nc)
  int             Nw;     // Number of poloidal grid-points on visualization grid (from Equilibrium.nc)
  double*         rf;     // Radial grid-points on visualization grid (from Equilibrium.nc)
  Array<double,2> RR;     // R coordinates of visualization grid-points (from Equilibrium.nc)
  Array<double,2> ZZ;     // Z coordinates of visualization grid-points (from Equilibrium.nc)
  Array<double,2> rvals;  // r values of visualization grid-points (from Equilibrium.nc)
  Array<double,2> thvals; // theta values of visualization grid-points (from Equilibrium.nc)

  Array<complex<double>,3> Psiuf;  // y components of Fourier transformed ideal eigenfunctions
  Array<complex<double>,3> Zuf;    // Z components of Fourier transformed ideal eigenfunctions
  Array<complex<double>,3> Psiuv;  // y components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,3> Zuv;    // Z components of unreconnected tearing eigenfunctions on visualization grid
   
  // --------------------
  // Vacuum solution data
  // --------------------
  double sa;                       // Edge magnetic shear

  Array<complex<double>,2> Pvac;   // Vacuum solution matrix
  Array<complex<double>,2> Qvac;   // Vacuum solution matrix
  Array<complex<double>,2> Rvac;   // Vacuum solution matrix
  Array<complex<double>,2> Svac;   // Vacuum solution matrix

  Array<complex<double>,2> Pdag;   // Hermitian conjugate of Pvac
  Array<complex<double>,2> Qdag;   // Hermitian conjugate of Qvac
  Array<complex<double>,2> Rdag;   // Hermitian conjugate of Rvac
  Array<complex<double>,2> Sdag;   // Hermitian conjugate of Svac

  Array<complex<double>,2> PRmat;  // Pdag * Rvac
  Array<complex<double>,2> PRher;  // Hermitian component of PRmat
  Array<complex<double>,2> PRant;  // Anti-Hermitian component of PRmat

  Array<complex<double>,2> QSmat;  // Qdag * Svac
  Array<complex<double>,2> QSher;  // Hermitian component of QSmat
  Array<complex<double>,2> QSant;  // Anti-Hermitian component of QSmat

  Array<complex<double>,2> PSmat;  // Pdag * Svac - Rdag * Qvac

  Array<complex<double>,2> QPmat;  // Qvac * Pdag
  Array<complex<double>,2> QPher;  // Hermitian component of QPmat
  Array<complex<double>,2> QPant;  // Anti-Hermitian component of QPmat

  Array<complex<double>,2> RSmat;  // Rvac * Sdag
  Array<complex<double>,2> RSher;  // Hermitian component of RSmat
  Array<complex<double>,2> RSant;  // Anti-Hermitian component of RSmat

  Array<complex<double>,2> SPmat;  // Pvac * Sdag - Qvac * Rdag

  Array<complex<double>,2> RPmat;  // RPmat * Rvac = Pvac
  Array<complex<double>,2> RPdag;  // Hermitian conjugate of RPmat
  Array<complex<double>,2> Hmat;   // No-wall vacuum response matrix: Hmat = (1/2) (RPmat + RPdag)

  Array<complex<double>,2> iRPmat; // iRPmat * Pvac = Rvac
  Array<complex<double>,2> iRPdag; // Hermitian conjugate of iRPmat
  Array<complex<double>,2> iHmat;  // Inverse no-wall vacuum response matrix: iHmat = (1/2) (iRPmat + iRPdag)
  
  // ------------------
  // Wall solution data
  // ------------------
  double                   bw;      // Relative wall radius (read from Equilibrium JSON file)
  double*                  rho;     // Wall scaling vector

  Array<complex<double>,2> Rwal;    // Wall solution matrix
  Array<complex<double>,2> Swal;    // Wall solution matrix

  Array<complex<double>,2> iImat;   // Rwal * iImat = Swal
  Array<complex<double>,2> iIher;   // Hermitian component of iImat
  Array<complex<double>,2> iIant;   // Anti-Hermitian component of iImat

  Array<complex<double>,2> PImat;   // PImat = Pvac * iImat - Qvac
  Array<complex<double>,2> RImat;   // RImat = Rvac * iImat - Svac
  Array<complex<double>,2> IRmat;   // IRmat = Rvac * iImat 

  Array<complex<double>,2> RPImat;  // RPImat * RImat = PImat
  Array<complex<double>,2> RPIdag;  // Hermitian conjugate of RPImat
  Array<complex<double>,2> Gmat;    // Perfect-wall vacuum response matrix: Gmat = (1/2) (RPImat + RPIdag)

  Array<complex<double>,2> iRPImat; // iRPImat * PImat = RImat
  Array<complex<double>,2> iRPIdag; // Hermitian conjugate of iRPImat
  Array<complex<double>,2> iGmat;   // Perfect-wall vacuum response matrix: iGmat = (1/2) (iRPImat + iRPIdag)

  Array<complex<double>,2> Bmat;    // Bmat * RImat = IRmat
  Array<complex<double>,2> Cmat;    // Cmat = (Gmat - Hmat) * Bmat
  Array<complex<double>,2> Cher;    // Hermitian component of Cmat
  Array<complex<double>,2> Cant;    // Anti-Hermitian component of Cmat

  // -----------------
  // ODE solution data
  // -----------------
  double*                  Rgrid; // Radial grid-points for diagnostics
  double*                  Pgrid; // PsiN values for grid-points for diagnostics
  Array<double,2>          Etest; // Energy test for solution vectors versus radius
  Array<double,2>          Pnorm; // Norms of y components of solution vectors versus radius
  Array<double,2>          Znorm; // Norms of Z components of solution vectors versus radius
  double*                  hode;  // Step-length versus radius
  double*                  eode;  // Truncation error versus radius
  Array<complex<double>,3> YYY;   // Solution vectors versus radius

  // -----------------
  // Ideal energy data
  // -----------------
  Array<complex<double>,3> Psii;    // y components of solutions launched from magnetic axis
  Array<complex<double>,3> Zi;      // Z components of solutions launched from magnetic axis
  Array<complex<double>,2> Wmat;    // Plasma ideal energy matrix
  Array<complex<double>,2> Wher;    // Hermitian component of Wmat
  Array<complex<double>,2> Want;    // Anti-Hermitian component of Wmat
  double*                  Wval;    // Eigenvalues of symmeterized W-matrix
  Array<complex<double>,2> Vmat;    // Vacuum ideal energy matrix
  Array<complex<double>,2> Vher;    // Hermitian component of Vmat
  Array<complex<double>,2> Vant;    // Anti-Hermitian component of Vmat
  double*                  Vval;    // Eigenvalues of symmeterized V-matrix
  Array<complex<double>,2> Umat;    // Total ideal energy matrix
  Array<complex<double>,2> Uher;    // Hermitian component of Umat
  Array<complex<double>,2> Uant;    // Anti-Hermitian component of Umat
  double*                  Uval;    // Eigenvalues of symmeterized U-matrix
  Array<complex<double>,2> Uvec;    // Eigenvectors of symmeterized U-matrix
  Array<complex<double>,2> Ures;    // Residuals of Uvec orthonormality matrix
  Array<complex<double>,3> Psie;    // y components of ideal eigenfunctions
  Array<complex<double>,3> Ze;      // Z components of ideal eigenfunctions
  double*                  deltaW;  // Total perturbed ideal potential energy 
  double*                  deltaWp; // Plasma contribution to perturbed ideal potential energy
  double*                  deltaWv; // Vacuum contribution to perturbed ideal potential energy
  Array<complex<double>,2> Psiy;    // y values on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> Xiy;     // Z values on plasma boundary associated with ideal eigenfunctions

  // ----
  // Misc
  // ----
   
 public:

  // ...............
  // In Vertical.cpp
  // ...............

  // Constructor
  Vertical ();
  // Destructor
  ~Vertical ();

  // Solve problem
  void Solve ();

 private:

  // ...............
  // In Vertical.cpp
  // ...............
  
  // Set mode numbers
  void SetModeNumbers ();
  // Deallocate memory
  void CleanUp ();

  // ......................
  // In ReadEquilibrium.cpp
  // ......................

  // Read equilibrium data from Outputs/Equilibrium/Equilibrium.nc
  void ReadEquilibrium ();
  // Calculate metric data at plasma boundary
  void CalculateMetricBoundary ();
  // Calculate metric data at wall
  void CalculateMetricWall ();

  // ..................
  // In Interpolate.cpp
  // ..................

  // Return value of PsiN
  double GetPsiN (double r);
  // Return value of f
  double Getf (double r);
  // Return value of p
  double Getp (double r);
  // Return value of pp
  double Getpp (double r);
  // Return value of ppp
  double Getppp (double r);
  // Return value of q
  double Getq (double r);
  // Return value of g2
  double Getg2 (double r);
  // Return value of g
  double Getg (double r);
  // Return value of s
  double Gets (double r);
  // Return value of s2
  double Gets2 (double r);
  // Return value of s0
  double Gets0 (double r);
  // Return value of S1
  double GetS1 (double r);
  // Return value of S2
  double GetS2 (double r);
  // Return value of S3
  double GetS3 (double r);
  // Return value of S4
  double GetS4 (double r);
  // Return value of S5
  double GetS5 (double r);
  // Return value of d Sig/dr
  double GetSigp (double r);
  // Return value of P1
  double GetP1 (double r);
  // Return value of P2
  double GetP2 (double r);
  // Return value of P3
  double GetP3 (double r);
  // Return value of Sig
  double GetSig (double r);
  // Return value of ne
  double Getne (double r);
  // Return value of Te
  double GetTe (double r);
  // Return value of nep
  double Getnep (double r);
  // Return value of Tep
  double GetTep (double r);
  // Return value of Hn
  double GetHn (int n, double r);
  // Return value of Hnp
  double GetHnp (int n, double r);
  // Return value of Vn
  double GetVn (int n, double r);
  // Return value of Vnp
  double GetVnp (int n, double r);

  // ...............
  // In Toroidal.cpp
  // ...............

  // Return associated Legendre function P^m_(n-1/2) (z)
  double ToroidalP (int m, int n, double z);
  // Return associated Legendre function Q^m_(n-1/2) (z)
  double ToroidalQ (int m, int n, double z);
  // Return derivative of associated Legendre function P^m_(n-1/2) (z)
  double ToroidaldPdz (int m, int n, double z);
  // Return derivative of associated Legendre function Q^m_(n-1/2) (z)
  double ToroidaldQdz (int m, int n, double z);
  // Return normalized associated Legendre function P^m_(n-1/2) (z)
  double NormToroidalP (int m, int n, double z);
  // Return normalized associated Legendre function Q^m_(n-1/2) (z)
  double NormToroidalQ (int m, int n, double z);
  // Return normalzied derivative of associated Legendre function P^m_(n-1/2) (z)
  double NormToroidaldPdz (int m, int n, double z);
  // Return normalized derivative of associated Legendre function Q^m_(n-1/2) (z)
  double NormToroidaldQdz (int m, int n, double z);
  // Return hyperbolic cosine of toroidal coordinate mu
  double GetCoshMu (double R, double Z);
  // Return toroidal coordinate eta
  double GetEta (double R, double Z);

  // .............
  // In Vacuum.cpp
  // .............

  // Calculate vacuum boundary matrices
  void GetVacuumBoundary ();
  // Calculate vacuum wall matrices
  void GetVacuumWall ();
  // Evaluate right-hand sides of vacuum odes
  void CashKarp45Rhs1 (double r, complex<double>* Y, complex<double>* dYdr) override;

  // ................
  // In Armadillo.cpp
  // ................

  // Solve linear system of equations A . X = B, for X, where all quantities are real rectangular matrices
  void SolveLinearSystem (Array<double,2> A, Array<double,2> X, Array<double,2> B);
  // Solve linear system of equations X . A = B, for X, where all quantities are real rectangular matrices
  void SolveLinearSystemTranspose (Array<double,2> A, Array<double,2> X, Array<double,2> B);
  // Solve linear system of equations A . X = B, for X, where all quantities are complex rectangular matrices
  void SolveLinearSystem (Array<complex<double>,2> A, Array<complex<double>,2> X, Array<complex<double>,2> B);
  // Solve linear system of equations X . A = B, for X, where all quantities are complex rectangular matrices
  void SolveLinearSystemTranspose (Array<complex<double>,2> A, Array<complex<double>,2> X, Array<complex<double>,2> B);
  // Solve linear system of equations A . X = B, for X, where A is a complex rectangular matrix, and X and B
  //  are complex vectors
  void SolveLinearSystem (Array<complex<double>,2> A, complex<double>* X, complex<double>* B);
  // Invert square complex matrix
  void InvertMatrix (Array<complex<double>,2> A, Array<complex<double>,2> invA);
  // Return eigenvalues and eigenvectors of Hermitian matix H
  void GetEigenvalues (Array<complex<double>,2> H, double* evals, Array<complex<double>,2> evecs);
  // Return eigenvalues of Hermitian matix H
  void GetEigenvalues (Array<complex<double>,2> H, double* evals);
  // Return eigenvalues and eigenvectors of pair of complex matrices
  void GetEigenvalues (Array<complex<double>,2> A, Array<complex<double>,2> B, complex<double>* evals, Array<complex<double>,2> evecs);
  // Find square root of Hermitian positive definite matrix
  void SquareRootMatrix (Array<complex<double>,2> A, Array<complex<double>,2> sqrtA);
  // Find inverse square root of Hermitian positive definite matrix
  void InvSquareRootMatrix (Array<complex<double>,2> A, Array<complex<double>,2> invsqrtA);

  // .............
  // In Matrix.cpp
  // .............

  // Get values of coupling matrices
  void GetMatrices (double r, 
		    Array<complex<double>,2> AAmmp, Array<complex<double>,2> BBmmp,
		    Array<complex<double>,2> NNmmp, Array<complex<double>,2> PPmmp);
  // Get values of coupling matrix elements
  void GetMatrices (double r, int m, int mp, 
		    complex<double>& Ammp, complex<double>& Bmmp,
		    complex<double>& Cmmp, complex<double>& Dmmp);

  // ...............
  // In ODESolve.cpp
  // ...............

  // Function to solve outer region odes
  void ODESolve ();
  // Launch solution vectors from magnetic axis
  void LaunchAxis (double r, Array<complex<double>,2> YY);
  // Integrate multiple solution vectors from given value of r to r = rx while performing nf fixups
  void SegmentFixup (double& r, double rx, int nf, Array<complex<double>,2> YY);
  // Integrate multiple solution vectors from given value of r to to r = rx
  void Segment (double& r, double rx, Array<complex<double>,2> YY);
  // Perform fixup of multiple solution vectors lauched from magnetic axis
  void Fixup (double r, Array<complex<double>,2> YY);
  // Evaluate right-hand sides of outer region odes
  void CashKarp45Rhs (double r, complex<double>* Y, complex<double>* dYdr) override;
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
  double EnergyTest (double r, Array<complex<double>,2> YY);
  // Calculate angular momentum flux associated with ith solution vector
  double GetEnergy (double r, Array<complex<double>,2> YY, int i);
  // Calculate norms of ith solution vector
  void GetNorms (Array<complex<double>,2> YY, int i, double& Pnorm, double &Znorm);

  // ............
  // In Ideal.cpp
  // ............

  // Calculate ideal stability
  void CalculateIdealStability ();

  // ................
  // In Visualize.cpp
  // ................
  
  // Calculate ideal eigenfunction visualization data
  void VisualizeEigenfunctions ();
 
  // .............
  // In Netcdf.cpp
  // .............

  // Read equilibrium data from netcdf file
  void ReadEquilibriumNetcdf ();
  // Write stability data to netcdf file
  void WriteNetcdf ();
};
