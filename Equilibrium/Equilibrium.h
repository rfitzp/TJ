// Equilibrium.h

// ###########################################################################################################

// Class to calculate inverse aspect-ratio expanded tokamak equilibrium.

// All lengths (except r) normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis toroidal magnetic field-strength).
// Radial coordinate, r, normalized to epsa * R_0, where epsa is inverse-aspect ratio of plasma.
// So r = 0 on magnetic axis and r = 1 at plasma/vacuum interface.

// Equilibrium magnetic flux-surfaces are defined parametrically as:

// R(r,w) = 1 - epsa r cosw + epsa^2 H1(r) + epsa^2 sum_{n=2,Ns} [Hn(r) cos(n-1)w + Vn(r) sin(n-1)w]
// Z(r,w) =     epsa r sinw                + epsa^2 sum_{n=2,Ns} [Hn(r) sin(n-1)w - Vn(r) cos(n-1)w]
//
// Here, R, phi, Z are cylindrical polar coordinates, r is a flux-surface label, and w is a poloidal angle.
// The class also uses the r, theta, phi (PEST) straight field-line coordinate system whose Jacobian is r R^2. 

// Edge shaping is specified as: Hna = Hn(1), Vna = Vn(1), etc.

// Inputs:
//  Inputs/Equilibrium.json - JSON file

// Outputs:
//  Outputs/Equilibrium/Equilibrium.nc

// Plotting scripts:
//  Plots/Equilibrium/*.py

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  GNU scientific library (https://www.gnu.org/software/gsl)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Author:
//  Richard Fitzpatrick,
//  Institute of Fusion Studies,
//  Department of Physics
//  University of Texas at Austin
//  rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/TJPaper/TJ.pdf

// ###########################################################################################################

#pragma once

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <netcdf>
#include "Utility.h"

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
    
// ############
// Class header
// ############
class Equilibrium : private Utility
{
 private:

  // -----------------
  // Calculation flags
  // -----------------
  int VIZ;  // Flag for calculating magnetic flux-surfaces (read from TJ JSON file)
  int SRC;  // Flag for reading lowest-order safety factor and pressure profiles from file:
            //  SRC = 0 - profiles are evaluated analytically:
            //
            //             q  = r^2 /f1
            //             f1 = [1 - (1 - r^2)^nu] /qc/nu
            //             p2 = pc (1 - r^2)^mu
            //
            //             qc is lowest-order safety-factor on magnetic axis.
            //             nu * qc is lowest-order safety-factor at plasma/vacuum interface.
            //
            //  SRC = 1 - profiles are read from Inputs/Profile.txt:
            //
            //              Format:
            //              
            //               NPTS
            //               r(0)      q(0)      p2(0)
            //               r(1)      q(1)      p2(1)
            //               ...       ...
            //               r(NPTS-1) q(NPTS-1) p2(NPTS-1)
            //
            //              Here, r is a cylindrical radial coordinate.
            //              Note that p2 is renormalized such that p2(0) = pc.
            //              qc is set to q(0), qa to q(NPTS-1), and nu to qa/qc.
            //
            //              LightEquilibrium is not called, so the true edge safety-factor
            //              will differ from qa
  
  // ------------------
  // Physics parameters
  // ------------------
  double         epsa;  // Inverse aspect-ratio of plasma (read from JSON file)
  double         qc;    // Lowest-order safety-factor on magnetic axis (read from JSON file)
  double         qa;    // Safety-factor at plasma/vacuum interface (read from JSON file)
  double         pc;    // Normalized plasma pressure on magnetic axis (read from JSON file)
  double         mu;    // Pressure peaking parameter (read from JSON file)
  vector<double> Hna;   // H2(1), H3(1), etc (read from JSON file)
  vector<double> Vna;   // V2(1), V3(1), etc (read from JSON file)
  double         nu;    // Toroidal current peaking parameter (determined from qa)

  double         bw;    // Relative wall radius (read from JSON file)

  int            ncoil; // Number of toroidal strands that make up RMP coils (strands assumed to be located at wall radius)
  double*        wcoil; // omega coordinates of RMP strands/M_PI (read from JSON file)
  double*        Icoil; // Toroidal currents flowing in RMP strands (read from JSON file)
  double*        Rcoil; // R coodinates of RMP strands 
  double*        Zcoil; // Z coodinates of RMP strands 

  // ------------------
  // Machine parameters
  // =-----------------
  double         R0;    // Major radius of magnetic axis (read from JSON file)
  double         B0;    // Toroidal magnetic field-strength on magnetic axis (read from JSON file)
  double         n0;    // Electron number density on magnetic axis (read from JSON file)
  double         alpha; // Electron number density profile parameter (read from JSON file)
  double         Teped; // Electron temperature at plasma boundnary (read from JSON file)
  double         neped; // Electron number at plasma boundnary (read from JSON file)

  // ------------------
  // Profile parameters
  // ------------------
  int                NPTS;        // Number of data points in profile file (Inputs/Profile.txt)
  double*            rin;         // Radial grid-points in profile file
  double*            qin;         // Lowest order safety-factor profile read from profile file
  double*            f1in;        // f1 profile derived from qin; f1in = r^2/qin
  double*            p2in;        // p2 profile read from profile file

  gsl_spline*        f1inspline;  // Interpolated f1in function
  gsl_spline*        p2inspline;  // Interpolated p2in function
  gsl_interp_accel*  f1inacc;     // Accelerator for interpolated f1in function
  gsl_interp_accel*  p2inacc;     // Accelerator for interpolated p2in function

  // ----------------------
  // Calculation parameters
  // ----------------------
  double eps; // Distance of closest approach to magnetic axis (read from JSON file)
  int    Ns;  // Number of shaping harmonics (read from JSON file)
  int    Nr;  // Number of radial grid-points for calculation purposes (read from JSON file)
  int    Nf;  // Number of radial grid-points for visualization purposes (read from JSON file)
  int    Nw;  // Number of angular grid-points for visualization purposes (read from JSON file)

  // ----------------
  // Calculation data
  // ----------------
  double*            rr;        // Radial grid-points
  double*            p2;        // Plasma pressure profile
  double*            f1;        // Lowest-order poloidal flux function
  double*            f3;        // Higher-order poloidal flux function
  double*            g2;        // Lowest-order toroidal flux function
  double*            q0;        // Lowest-order safety-factor
  double*            q2;        // Higher-order safety-factor
  double*            It;        // Toroidal plasma current
  double*            Ip;        // Poloidal plasma current
  double*            Jt;        // Radial derivative of toroidal plasma current
  double*            Jp;        // Radial derivative of poloidal plasma current
  double*            pp;        // Radial derivative of plasma pressure
  double*            ppp;       // Second radial derivative of plasma pressure 
  double*            qq;        // First radial derivative of safety-factor times r
  double*            qqq;       // Radial derivative of qq times r
  double*            s;         // Higher-order magnetic shear:  s  = r q2'/q2
  double*            s2;        // Derivative of magnetic shear: s2 = r^2 q2''/q2
  double*            s0;        // Lowest-order magnetic shear (s0=2 at boundary)
  double*            S1;        // First shaping function
  double*            S2;        // Second shaping function
  double*            S3;        // Third shaping function
  double*            S4;        // Fourth shaping function
  double*            P1;        // First profile function: (2-s) /q2
  double*            P1a;       // Lowest-order first profile function: (2-s0) /q0 
  double*            P2;        // Second profile function: r dP1/dr
  double*            P2a;       // Lowest-order second profile function: r dP1a/dr
  double*            P3;        // Third profile function
  double*            P3a;       // Auxillary third profile function
  double*            ff;        // f profile
  double*            ggr2;      // <|nabla r|^2> profile
  double*            RR2;       // <R^2> profile
  double*            IR2;       // <|nabla r|^2/R^2> profile
  double*            Psi;       // Psi(r) array
  double*            PsiN;      // PsiN(r) array
  double*            Tf;        // Toroidal flux-function
  double*            mu0P;      // mu_0 times pressure
  double*            DI;        // Ideal Mercier index
  double*            DR;        // Resistive Mercier index
  double*            Te;        // Electron temperature
  double*            ne;        // Electron number density  
  double*            Tep;       // Radial derivative of electron temperature
  double*            nep;       // Radial derivative of electron number density
   
  Array<double,2>    HHfunc;    // Horizontal shaping functions
  Array<double,2>    VVfunc;    // Vertical shaping functions
  Array<double,2>    HPfunc;    // Radial derivatives of horizontal shaping functions
  Array<double,2>    VPfunc;    // Radial derivatives of vertical shaping functions
  double*            Lfunc;     // Relabelling function
 
  double             amean;     // Mean minor radius
  double             li;        // Normalized plasma self-inductance
  double             betat;     // Plasma toroidal beta
  double             betap;     // Plasma poloidal beta
  double             betaN;     // Plasma normal beta
  
  gsl_spline*        Itspline;  // Interpolated It function
  gsl_spline*        Ipspline;  // Interpolated Ip function
  gsl_spline*        g2spline;  // Interpolated g2(r) function
  gsl_spline**       HHspline;  // Interpolated horizontal shaping functions versus r
  gsl_spline**       VVspline;  // Interpolated vertical shaping functions versus r
  gsl_spline**       HPspline;  // Interpolated radial derivatives of horizontal shaping functions versus r
  gsl_spline**       VPspline;  // Interpolated radial derivatives of vertical shaping functions versus r
  gsl_spline*        Lspline;   // Interpolated relabelling function versus r
  gsl_spline*        wspline;   // Interpolated omega(theta) function
  gsl_spline*        Rspline;   // Interpolated R(theta) function
  gsl_spline*        Zspline;   // Interpolated Z(theta) function
  gsl_spline*        fspline;   // Interpolated f function versus r
  gsl_spline*        q2spline;  // Interpolated q2 versus r
  gsl_spline*        gr2spline; // Interpolated <|nabla r|^2> function versus r
  gsl_spline*        R2spline;  // Interpolated <R^2> function versus r
  gsl_spline*        I2spline;  // Interpolated <|nabla r|^2/R^2> function versus r
  gsl_spline*        sspline;   // Interpolated s function versus r
  gsl_spline*        qspline;   // Interpolated q function versus r

  // ------------------
  // Vizualization data
  // ------------------
  gsl_interp_accel*  Itacc;     // Accelerator for interpolated It function
  gsl_interp_accel*  Ipacc;     // Accelerator for interpolated Ip function
  gsl_interp_accel*  g2acc;     // Accelerator for interpolated g2 function
  gsl_interp_accel** HHacc;     // Accelerators for interpolated horizontal shaping functions
  gsl_interp_accel** VVacc;     // Accelerators for interpolated vertical shaping functions
  gsl_interp_accel** HPacc;     // Accelerators for interpolated radial derivatives of horizontal shaping functions
  gsl_interp_accel** VPacc;     // Accelerators for interpolated radial derivatives of vertical shaping functions
  gsl_interp_accel*  Lacc;      // Accelerator for interpolated relabelling function
  gsl_interp_accel*  wacc;      // Accelerator for interpolated omega function
  gsl_interp_accel*  Racc;      // Accelerator for interpolated R(theta)
  gsl_interp_accel*  Zacc;      // Accelerator for interpolated Z(theta)
  gsl_interp_accel*  facc;      // Accelerator for interpolated f function
  gsl_interp_accel*  q2acc;     // Accelerator for interpolated q2 function
  gsl_interp_accel*  gr2acc;    // Accelerator for interpolated <|nabla r|^2> function
  gsl_interp_accel*  R2acc;     // Accelerator for interpolated <R^2> function
  gsl_interp_accel*  I2acc;     // Accelerator for interpolated <R^2> function
  gsl_interp_accel*  sacc;      // Accelerator for interpolated s function
  gsl_interp_accel*  qacc;      // Accelerator for interpolated q function

  Array<double,2>    RR;        // R coodinates of magnetic flux-surfaces for visualization purposes (uniform theta grid)
  Array<double,2>    ZZ;        // Z coodinates of magnetic flux-surfaces for visualization purposes (uniform theta grid)
  Array<double,2>    dRdr;      // dR/dr values on visualization grid
  Array<double,2>    dRdt;      // (dR/dtheta)/r values on visualization grid
  Array<double,2>    dZdr;      // dZ/dr values on visualization grid
  Array<double,2>    dZdt;      // (dZ/dtheta)/r values on visualization grid
  Array<double,2>    Jac;       // Jacobian
  Array<double,2>    Jax;       // Analytic Jacobian
  Array<double,2>    RRw;       // R coodinates of magnetic flux-surfaces for visualization purposes (uniform omega grid)
  Array<double,2>    ZZw;       // Z coodinates of magnetic flux-surfaces for visualization purposes (uniform omega grid)
  Array<double,2>    rvals;     // r values on magnetic flux-surfaces for visualization purposes
  Array<double,2>    thvals;    // theta values on magnetic flux-surfaces for visualization purposes
  Array<double,2>    wvals;     // omega values on magnetic flux-surfaces for visualization purposes

  // -------------
  // Boundary data
  // -------------
  double*            Rbound;    // R values on plasma boundary
  double*            Zbound;    // Z values on plasma boundary
  double*            tbound;    // theta values on plasma boundary
  double*            wbound;    // omega values on plasma boundary
  double*            tbound0;   // Preliminary theta values on plasma boundary
  double*            wbound0;   // Preliminary omega values on plasma boundary
  double*            dRdtheta;  // dR/dtheta values on plasma boundary
  double*            dZdtheta;  // dZ/dtheta values on plasma boundary
  double*            R2b;       // R^2 values on plasma boundary
  double*            grr2b;     // |nabla r|^2 values on plasma boundary

  // ---------
  // Wall data
  // ---------
  double*            Rwall;     // R values on wall
  double*            Zwall;     // Z values on wall
  double*            wwall;     // omega values on wall

  // -------------------------
  // Synthetic diagnostic data
  // -------------------------
  double  tilt;   // tilt of central chord (degrees) (read from JSON file)
  double* req;    // r values on central chord
  double* weq;    // omega values on central chord
  double* teq;    // theta values on central chord
  double* Req;    // R values on central chord
  double* Zeq;    // Z values on central chord
  double* BReq;   // B_parallel values on central chord
  double* neeq;   // n_e values on central chord
  double* Teeq;   // T_e values on central chord
  double* dRdreq; // dRdr values on central chord
  double* dRdteq; // dRdt values on central chord
  double* dZdreq; // dZdr values on central chord
  double* dZdteq; // dZdt values on central chord
  double* Weeq;   // First harmonic cyclotron frequency on central chord
  double* wLeq;   // Lower cut-off frequency on central chord
  double* wUeq;   // Upper cut-off frequency on central chord
  double* wUHeq;  // Upper-hybrid frequency on central chord

  // ---------------
  // EFIT parameters
  // ---------------
  int                EFIT;       // Flag for calculating EFIT data (read from JSON file)
  int                NRBOX;      // Number of R grid-points (read from JSON file)
  int                NZBOX;      // Number of Z grid-points (read from JSON file)
  double             rc;         // Flux-surfaces calculated accurately up to r = rc (read from JSON file)
  int                NPBOUND;    // Number of boundary points (= NW+1)
  int                NLIMITER;   // Number of limiter points (= 5)
  double             RBOXLFT;    // Left-hand coordinate of R box (deduced from boundary values)
  double             RBOXLEN;    // Length of R box (deduced from boundary values)
  double             ZOFF;       // Offset of centroid of Z box (deduced from boundary values)
  double             ZBOXLEN;    // Length of Z box (deduced from boundary values)
  double             R0EXP;      // Major radius of magnetic axis (read from JSON file)
  double             B0EXP;      // Toroidal magnetic field-strength on magnetic axis (read from JSON file)
  double             RAXIS;      // R coordinate of magnetic axis (= R0)
  double             ZAXIS;      // Z coordinate of magnetic axis (= Z0)
  double             PSIAXIS;    // PSI value on magnetic axis
  double             PSIBOUND;   // PSI value on plasma boundary (= 0)
  double             CURRENT;    // Toroidal plasma current
  double*            PSI;        // Equally-spaced PSI array
  double*            PSIN;       // PSI_N array
  double*            rPSI;       // r values that coincide with PSI values
  double*            PSIr;       // PSI values that coincide with r values
  double*            T;          // Toroidal magnetic flux evaluated on PSI grid
  double*            TTp;        // T T' evaluated on PSI grid
  double*            P;          // P evaluated on PSI grid
  double*            Pp;         // P' evaluated on PSI grid
  double*            Q;          // Safety-factor evaluated on PSI grid
  double*            RBOUND;     // R values on plasma boundary
  double*            ZBOUND;     // Z values on plasma boundary
  double*            RLIMITER;   // R values on limiter
  double*            ZLIMITER;   // Z values on limiter
  double*            RGRID;      // R grid-points
  double*            ZGRID;      // Z grid-points
  double*            PSIRZ;      // PSI evaluated on R, Z grid
  double*            rRZ;        // r evaluated on R, Z grid
  double*            wRZ;        // w evaluated on R, Z grid
  double*            cwRZ;       // cos(w) evaluated on R, Z grid
  double*            swRZ;       // sin(w) evaluated on R, Z grid

  gsl_spline*        rPsispline; // Interpolated r function versus Psi
  gsl_spline*        PSIrspline; // Interpolated PSI function versus r
  gsl_interp_accel*  rPsiacc;    // Accelerator for interpolated r function versus Psi
  gsl_interp_accel*  PSIracc;    // Accelerator for interpolated PSI function versus r

  // ----
  // Misc
  // ----
  int rhs_chooser;

 public:

  // ..................
  // in Equilibrium.cpp
  // ..................
  
  // Constructor
  Equilibrium ();

  // Destructor
  ~Equilibrium ();

  // Find nu value that corresponds to edge safety-factor value read from JSON file
  void Setnu ();

  // Solve problem
  void Solve ();

private:

  // .............
  // in Netcdf.cpp
  // .............

  // Write equilibrium data to netcdf file
  void WriteNetcdf (double sa);

  // ...........
  // in EFIT.cpp
  // ...........
  
  // Calculate EFIT data
  void CalculateEFIT ();

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
  // Return q
  double Getq (double r);
  // Return s
  double Gets (double r);
  // Return DI
  double GetDI (double r);
  // Return DR
  double GetDR (double r);
  // Return Te
  double GetTe (double r);
  // Return ne
  double Getne (double r);
  // Return Tep
  double GetTep (double r);
  // Return nep
  double Getnep (double r);
  // Return Hnp
  double GetHnp (int n, double r);

  // Return relabelling parameter
  double GetL (double r, int order);
  // Return R
  double GetR (double r, double w, int order);
  // Return dRdr
  double GetdRdr (double r, double w, int order);
  // Return dRdw
  double GetdRdw (double r, double w, int order);
  // Return Z
  double GetZ (double r, double w, int order);
  // Return dZdr
  double GetdZdr (double r, double w, int order);
  // Return dZdw
  double GetdZdw (double r, double w, int order);
  // Return theta (r, omega) function
  double Gettheta (double r, double w, int order);
  // Return dtheta (r, omega)/domega function
  double Getdthetadomega (double r, double w, int order);
  // Return R2
  double GetR2 (double r, double t);
  // Return H_j in vacuum
  double GetHHvac (int j, double r);
  // Return H_j' in vacuum
  double GetHPvac (int j, double r);
  // Return V_j in vacuum
  double GetVVvac (int j, double r);
  // Return V_j' in vacuum
  double GetVPvac (int j, double r);
  // Return Psi in vacuum
  double GetPSIvac (double r);
  // Return f_R
  double Getf_R (double r, double w);
  // Return f_Z
  double Getf_Z (double r, double w);
  // Return |nabla r|^2
  double Getgrr2 (double r, double t);
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double* y, double* dydx) override;
};
