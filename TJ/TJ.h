// TJ.h

// ###################################################################################

// Class to calculate tearing stability matrix and tearing eigenfunctions in an
// inverse aspect-ratio expanded tokamak equilibrium.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field).
// Radial coordinate, r, normalized to epsa * R_0, where epsa is inverse-aspect ratio.
// So r = 0 corresponds to magnetic axis, and r = 1 to plasma/vacuum interface.

// See Equilibrium.h for description of equilibrium.
// Program assumes monotonic safety-factor profile.

// Inputs:
//  Inputs/TJ.json - JSON file

// Outputs:
//  Outputs/TJ/TJ,nc

// Plotting scripts:
//  Plots/TJ/*.py

// Class uses following external libraries:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)
//  Armadillo library      (https://arma.sourceforge.net)

// Author:
// Richard Fitzpatrick,
// Institute of Fusion Studies,
// Department of Physics,
// University of Texas at Austin,
// rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ/

// Documentation: ../Documentation/TJPaper/TJ.pdf

// ###################################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <nlohmann/json.hpp>
#include <netcdf>
#include <armadillo>

using namespace blitz;
using           json = nlohmann::json;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace arma;

// ############
// Class header
// ############
class TJ
{
 private:

  // ----------------------
  // Calculation parameters
  // ----------------------
  int    NTOR;    // Toroidal mode number (read from JSON file)
  int    MMIN;    // Minimum poloidal mode number included in calculation (read from JSON file)
  int    MMAX;    // Maximum poloidal mode number included in calculation (read from JSON file)

  double EPS;     // Solutions launched from magnetic axis at r = EPS (read from JSON file)
  double DEL;     // Distance of closest approach to rational surface is DEL (read from JSON file)
  int    NFIX;    // Number of fixups performed (read from JSON file)
  int    NDIAG;   // Number of radial grid-points for diagnostics (read from JSON file)
  double NULC;    // Use zero pressure jump conditions when |nu_L| < NULC (read from JSON file)
  int    ITERMAX; // Maximum number of iterations used to determine quantities at rational surface (read from JSON file)
  int    FREE;    // Flag for free/fixed boundary calculation (read from JSON file)

  double EPSF;    // Step-length for finite difference determination of derivatives

  // ----------------------------
  // Layer calculation parameters
  // ----------------------------
  double B0;      // On-axis toroidal magnetic field-strength (T) (read from JSON file)
  double R0;      // Plasma minor radius (m) (read from JSON file)
  double n0;      // On-axis electron number density (m^-3) (read from JSON file)
  double alpha;   // Assumed electron number density profile: n0 (1 - r^2)^alpha (read from JSON file)
  double Zeff;    // Effective ion charge number (read from JSON file)
  double Mion;    // Ion charge number (read from JSON file)
  double Chip;    // Perpendicular momentum/energy diffusivity (m^2/s) (read from JSON file)
  double Teped;   // Electron temperature at edge of plasma (eV) (read from JSON file)
  double apol;    // Plasma minor radius (m)

  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;     // Integration accuracy (read from JSON file)
  double h0;      // Initial step-length (read from JSON file)
  double hmin;    // Minimum step-length (read from JSON file)
  double hmax;    // Maximum step-length (read from JSON file)
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Integration error calculation flag
  
  // ---------------------------------------------------------------
  // Equilibrium data (read from Outputs/Equilibrium/Equilibrium.nc)
  // ---------------------------------------------------------------
  double             epsa;      // Inverse aspect-ratio of plasma
  int                Ns;        // Number of shaping harmonics
  int                Nr;        // Number of radial grid-points

  double*            rr;        // Radial grid-points
  double*            g2;        // Second-order toroidal flux
  double*            p2;        // Second-order plasma pressure
  double*            pp;        // First radial derivative of second-order plasma pressure
  double*            ppp;       // Second radial derivative of second-order plasma pressure 
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

  gsl_spline*        g2spline;  // Interpolated p2 function
  gsl_spline*        p2spline;  // Interpolated p2 function
  gsl_spline*        ppspline;  // Interpolated pp function
  gsl_spline*        pppspline; // Interpolated ppp function
  gsl_spline*        qspline;   // Interpolated q function
  gsl_spline*        sspline;   // Interpolated s function
  gsl_spline*        s2spline;  // Interpolated s2 function
  gsl_spline*        S1spline;  // Interpolated S1 function
  gsl_spline*        P1spline;  // Interpolated P1 function
  gsl_spline*        P2spline;  // Interpolated P2 function
  gsl_spline*        P3spline;  // Interpolated P3 function

  gsl_interp_accel*  g2acc;     // Accelerator for interpolated g2 function
  gsl_interp_accel*  p2acc;     // Accelerator for interpolated p2 function
  gsl_interp_accel*  ppacc;     // Accelerator for interpolated pp function
  gsl_interp_accel*  pppacc;    // Accelerator for interpolated ppp function
  gsl_interp_accel*  qacc;      // Accelerator for interpolated q function
  gsl_interp_accel*  sacc;      // Accelerator for interpolated s function
  gsl_interp_accel*  s2acc;     // Accelerator for interpolated s2 function
  gsl_interp_accel*  S1acc;     // Accelerator for interpolated S1 function
  gsl_interp_accel*  P1acc;     // Accelerator for interpolated P1 function
  gsl_interp_accel*  P2acc;     // Accelerator for interpolated P2 function
  gsl_interp_accel*  P3acc;     // Accelerator for interpolated P3 function

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

  gsl_interp_accel*  Rrzacc;    // Accelerator for interpolated R2grgz function
  gsl_interp_accel*  Rreacc;    // Accelerator for interpolated R2grge function
  gsl_interp_accel*  Rbacc;     // Accelerator for interpolated R function
  gsl_interp_accel*  Zbacc;     // Accelerator for interpolated Z function

  // --------------------
  // Plasma boundary data
  // --------------------
  double* tbound; // theta values on plasma boundary
  double* Rbound; // R values on plasma boundary
  double* Zbound; // Z values on plasma boundary
  double* dRdthe; // dR/dtheta values on plasma boundary
  double* dZdthe; // dZ/dtheta values on plasma boundary

  // --------------------
  // Vacuum solution data
  // --------------------
  double                   sa;   // Edge magnetic shear
  Array<complex<double>,2> Pvac; // Vacuum solution matrix
  Array<complex<double>,2> Pdag; // Hermitian conjugate of Pvac
  Array<complex<double>,2> Rvac; // Vacuum solution matrix
  Array<complex<double>,2> Amat; // Pdag * Rvac
  Array<complex<double>,2> Aher; // Hermitian component of Amat
  Array<complex<double>,2> Aant; // Anti-Hermitian component of Amat
  Array<complex<double>,2> Rmat; // Pdag * Rmat = Aher
  Array<complex<double>,2> Rdag; // Hermitian conjugate of Rmat
  Array<complex<double>,2> Hmat; // Vacuum response matrix: Rdag * Hmat = Pdag

  // ---------------------
  // Rational surface data
  // ---------------------
  int     nres;    // Number of rational magnetic flux-surfaces
  int*    mres;    // Poloidal mode numbers at rational surfaces
  double* qres;    // Safety-factors at rational surfaces
  double* rres;    // Minor radii of rational surfaces
  double* qerr;    // Residual in determination of rational surface
  double* sres;    // Magnetic shears at rational surfaces
  double* DIres;   // DI values at rational surfaces
  double* DRres;   // DR values at rational surfaces
  double* nuLres;  // Mercier indices of large solution at rational surfaces
  double* nuSres;  // Mercies indices of small solution at rational surfaces
  int*    Jres;    // Index of resonant poloidal harmonic at rational surfaces

  // -------------------
  // Resonant layer data
  // -------------------
  double* Teres;  // Electron temperature (eV)
  double* S13res; // Cube root of Lundquist number
  double* taures; // Magnetic reconnection timescale S^(1/3) tauH (s)
  double* ieres;  // Minus ratio of electron diamagnetic frequency to ion diamagnetic frequency
  double* QEres;  // Normalized ExB frequency
  double* Qeres;  // Normalized electron diamagnetic frequency
  double* Qires;  // Normalized ion diamagnetic frequency
  double* Dres;   // Normalized ion sound radius
  double* Pmres;  // Magnetic Prandtl number for perpendicular momentum diffusion
  double* Peres;  // Magnetic Prandtl number for perpendicular energy diffusion
  double* Dcres;  // Critical Delta' for instability

  // ----------------
  // Mode number data
  // ----------------
  double   ntor; // Toroidal mode number
  int      J;    // Number of poloidal harmonics included in calculation
  int      K;    // Number of solution vectors: K = J + nres
  int*     MPOL; // Poloidal mode numbers of included poloidal harmonics
  double*  mpol; // Poloidal mode numbers of included poloidal harmonics

  // -----------------
  // ODE Solution data
  // -----------------
  double*                  Rgrid; // Radial grid-points for diagnostics
  Array<double,2>          Ttest; // Torque test for solution vectors versus radius
  Array<double,2>          Pnorm; // Norms of Psi components of solution vectors versus radius
  Array<double,2>          Znorm; // Norms of Z components of solution vectors versus radius
  double*                  hode;  // Step-length versus radius
  double*                  eode;  // Truncation error versus radius
  Array<complex<double>,3> YYY;   // Solution vectors versus radius
  Array<complex<double>,2> Pi;    // Reconnected fluxes at rational surfaces associated with solution vectors
  Array<complex<double>,1> dPi;   // Current sheets at rational surfaces associated with solution vectors

  // -------------------------------------
  // Tearing-mode dispersion relation data
  // -------------------------------------
  Array<complex<double>,2> Psia;  // Values of Psi at plasma boundary due to continuous solutions launched from magnetic axis
  Array<complex<double>,2> Za;    // Values of Z at plasma boundary due to continuous solutions launched from magnetic axis
  Array<complex<double>,2> Pia;   // Values of reconnected fluxes at rational surfaces due to continuous solutions
                                  //  launched from magnetic axis
  Array<complex<double>,2> Psis;  // Values of Psi at plasma boundary due to small solutions launched from rational surfaces
  Array<complex<double>,2> Zs;    // Values of Z at plasma boundary due to small solutions launched from rational surfaces
  Array<complex<double>,2> Pis;   // Values of reconnected fluxes at rational surfaces due to small solutions
                                  //  launched from rational surfaces
  Array<complex<double>,2> Xmat;  // X-matrix
  Array<complex<double>,2> Ymat;  // Y-matrix
  Array<complex<double>,2> Omat;  // Omega-matrix
  Array<complex<double>,2> Fmat;  // Inductance matrix
  Array<complex<double>,2> Fher;  // Symmeterized F-matrix
  double*                  Fval;  // Eigenvalues of symmeterized F-matrix
  Array<complex<double>,2> Fvec;  // Eigenvectors of symmeterized F-matrix
  Array<complex<double>,2> Emat;  // Tearing stability matrix
  Array<complex<double>,3> Psif;  // Psi components of fully reconnected tearing eigenfunctions
  Array<complex<double>,3> Zf;    // Z componnents of fully reconnected tearing eigenfunctions
  Array<complex<double>,3> Psiu;  // Psi components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> Zu;    // Z components of unreconnected tearing eigenfunctions
  Array<double,2>          Tf;    // Torques associated with fully reconnected eigenfunctions
  Array<double,2>          Tu;    // Torques associated with unreconnected eigenfunctions
  Array<double,3>          Tfull; // Torques associated with pairs of fully reconnected eigenfunctions
  Array<double,3>          Tunrc; // Torques associated with pairs of unreconnected eigenfunctions

  // -----------------------------------
  // Resonant magnetic perturbation data
  // -----------------------------------
  int                      ncoil;   // Number of toroidal strands that make up RMP coils 
  double*                  Rcoil;   // R coodinates of strands (read from JSON file)
  double*                  Zcoil;   // Z coodinates of strands (read from JSON file)
  double*                  Icoil;   // Toroidal currents flowing in strands (read from JSON file)
  complex<double>*         Psix;    // RMP perturbation at plasma boundary
  complex<double>*         Xi;      // RMP response vector
  complex<double>*         Upsilon; // RMP response vector
  complex<double>*         Lambda;  // RMP response vector
  complex<double>*         Chi;     // RMP drive at rational surfaces
  Array<complex<double>,2> Psirmp;  // Psi component of ideal RMP response eigenfunction
  Array<complex<double>,2> Zrmp;    // Z component of ideal RMP response eigenfunction
  complex<double>*         Psixs;   // Psix on plasma boundary
  complex<double>*         Psirmps; // Psirmp on plasma boundary

  // --------------------
  // Ideal stability data
  // --------------------
  int                      XiFlag;  // Flag for using Xi, rather than Psi, as ideal eigenfunction basis (read from JSON file)
  Array<complex<double>,3> Psii;    // Psi components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> Zi;      // Z components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> Xii;     // Xi components of ideal solutions launched from magnetic axis
  Array<complex<double>,2> Ji;      // Poloidal harmonics of current on plasma boundary associated with
                                    //  ideal solutions launched from magnetic axis
  Array<complex<double>,2> Wmat;    // Plasma ideal energy matrix
  Array<complex<double>,2> Vmat;    // Vacuum ideal energy matrix
  Array<complex<double>,2> Umat;    // Total ideal energy matrix
  Array<complex<double>,2> Uher;    // Hermitian component of Umat
  Array<complex<double>,2> Uant;    // Anti-Hermitian component of Umat
  double*                  Uval;    // Eigenvalues of symmeterized U-matrix
  Array<complex<double>,2> Uvec;    // Eigenvectors of symmeterized U-matrix
  Array<complex<double>,2> Ures;    // Residuals of Uvec orthonormaility matrix
  Array<complex<double>,3> Psie;    // Psi components of ideal eigenfunctions
  Array<complex<double>,3> Ze;      // Z components of ideal eigenfunctions
  Array<complex<double>,3> Xie;     // Xi components of ideal eigenfunctions
  Array<complex<double>,2> Je;      // Poloidal harmomics of current on plasma boundary associated with ideal eigenfunctions
  double*                  deltaW;  // Total perturbed ideal potential energy 
  double*                  deltaWp; // Plasma contribution to perturbed ideal potential energy
  double*                  deltaWv; // Plasma contribution to perturbed ideal potential energy
  Array<complex<double>,2> Psiy;    // Psi values on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> Jy;      // Current on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> Xiy;     // Xi values on plasma boundary associated with ideal eigenfunctions
  complex<double>*         gammax;  // Expansion of Psi_x at boundary in ideal eigenfunctions
  complex<double>*         gamma;   // Expansion of Psi_rmp at boundary in ideal eigenfunctions

  // ------------------------------------------------
  // Visualization of tearing eigenfunctions and RMPs
  // ------------------------------------------------
  int                      Nf;     // Number of radial grid-points on visualization grid
  int                      Nw;     // Number of angular grid-points on visualization grid
  double*                  rf;     // Radial grid-points on visualization grif
  Array<double,2>          RR;     // R coordinates of visualization grid-points
  Array<double,2>          ZZ;     // Z coordinates of visualization grid-points
  Array<double,2>          rvals;  // r values of visulalization grid-points
  Array<double,2>          thvals; // theta values of visualization grid-points
  Array<complex<double>,3> Psiuf;  // Psi components of Fourier-transformed unreconnected tearing eigenfunctions 
  Array<complex<double>,3> Zuf;    // Z components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> Psiuv;  // Psi components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,3> Zuv;    // Z components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,2> Psirf;  // Psi components of Fourier-transformed ideal RMP response eigenfunction
  Array<complex<double>,2> Zrf;    // Z components of Fourier-transformed ideal RMP response eigenfunction
  Array<complex<double>,2> Psirv;  // Psi components of ideal RMP response eigenfunction
  Array<complex<double>,2> Zrv;    // Z components of ideal RMP response eigenfunction

  // -----------------------
  // Root finding parameters
  // -----------------------
  double Eta;     // Minimum magnitude of f at root f(x) = 0
  int    Maxiter; // Maximum number of iterations

  // ----------------------------
  // Cash-Karp RK4/RK5 parameters
  // ----------------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int    count, rhs_chooser;
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
  // Deallocate memory
  void CleanUp ();
  // Read JSON file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (char* filename);
  // Open file for reading
  FILE* OpenFiler (char* filename);

  // ......................
  // In ReadEquilibrium.cpp
  // ......................

  // Read equilibrium data from Inputs/Equilibrium.nc
  void ReadEquilibrium ();
  // Read RMP coil data
  void ReadCoils ();
  // Calculate metric data at plasma boundary
  void CalculateMetric ();
    
  // .............
  // In Vacuum.cpp
  // .............

  // Calculate vacuum matrices
  void GetVacuum ();
  // Evaluate right-hand sides of vacuum odes
  void Rhs1 (double r, complex<double>* Y, complex<double>* dYdr);
  // Evaluate eta for RMP coil calculation
  double Geteta (double R, double Z, double Rp, double Zp);
  // Evaluate G for RMP coil calculation
  double GetG (double R, double Z, double Rp, double Zp);
  
  // ...............
  // In Rational.cpp
  // ...............

  // Target function for finding rational surfaces
  double Feval (double r);
  // Find rational surfaces
  void FindRational ();
  // Get resonant layer data
  void GetLayerData ();

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
  // Calculate resonant magnetic perturbation data
  void CalculateResonantMagneticPerturbation ();
  // Calculate unreconnected eigenfunction and RMP response visualization data
  void VisualizeEigenfunctions ();
  // Calculate ideal stability
  void CalculateIdealStability ();
   
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
  // Return value of DR
  double GetDR (double r);

  // ...............
  // In ZeroFind.cpp
  // ...............

  // Routine to find approximate root of F(x) = 0 using Ridder's method
  double RootFind ();
  // Ridder's method for finding root of F(x) = 0
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
  // Advance set of coupled first-order o.d.e.s by single step using adaptive
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Adaptive1 (int neqns, double& x, complex<double>* y, double& h, 
			    double& t_err, double acc, double S, double T, int& rept,
			    int maxrept, double h_min, double h_max, int flag, 
			    int diag, FILE* file);
  // Advance set of coupled first-order o.d.e.s by single step using fixed
  //  step-length Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Fixed1 (int neqns, double& x, complex<double>* y, complex<double>* err, double h);
 
  // ................
  // In Armadillo.cpp
  // ...............

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
  // Return hyperbolic cosine of toroidal coordinate mu
  double GetCoshMu (double R, double Z);
  // Return toroidal coordinate eta
  double GetEta (double R, double Z);

  // .............
  // In Netcdf.cpp
  // .............

  // Read equilibrium data from netcdf file
  void ReadNetcdf ();
  // Write stability data to netcdf file
  void WriteNetcdf ();
};
