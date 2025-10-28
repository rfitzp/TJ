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
//  Inputs/Equilibrium.json  - JSON file
//  Inputs/TJ.json           - JSON file
//  Inputs/Island.json       - JSON file

// Outputs:
//  Outputs/TJ/TJ.nc

// Plotting scripts:
//  Plots/TJ/*.py

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

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/

// ###################################################################################

#pragma once

#include "Utility.h"
#include "Equilibrium.h"
#include "Layer.h"
#include "Island.h"
#include "ECE.h"

// ############
// Class header
// ############
class TJ : private Utility
{
 private:

  // -----------------
  // Calculation flags
  // -----------------
  int SRC;     // Flag for reading profile data from file (read from Equilibrium JSON file)
  int EQLB;    // Flag for equilibrium calculation only (read from TJ JSON file)
  int FREE;    // Flag for free/fixed boundary calculation (read from TJ JSON file)
               //  FREE > 0 - perform no-wall tearing mode calculation
               //  FREE = 0 - perform perfect-wall tearing mode calculation
               //  FREE < 0 - perform fixed-boundary tearing mode calculation
  int FVAL;    // Flag for calculating eigenvalues and eigenvectors of F-matrix (read from TJ JSON file)
  int RMP;     // Flag for resonant magnetic perturbation calculation (read from TJ JSON file)
  int VIZ;     // Flag for generating unreconnected eigenfunction visualization data (read from TJ JSON file)
  int IDEAL;   // Flag for ideal calculation (read from TJ JSON file)
  int XI;      // Flag for using Xi, rather than Psi, as ideal eigenfunction basis (read from TJ JSON file)
  int LAYER;   // Flag for layer calculation (read from TJ JSON file)
  int TEMP;    // Flag for perturbed temperature calculation (read from TJ JSON file)

  // ----------------------
  // Calculation parameters
  // ----------------------
  int            NTOR;    // Toroidal mode number (read from TJ JSON file)
  int            MMIN;    // Minimum poloidal mode number included in calculation (read from TJ JSON file)
  int            MMAX;    // Maximum poloidal mode number included in calculation (read from TJ JSON file)
  vector<double> ISLAND;  // Island widths/displacements (divided by minor radius) used to regularize perturbed magnetic field (read from TJ JSON file)

  double         EPS;     // Solutions launched from magnetic axis at r = EPS (read from TJ JSON file)
  double         DEL;     // Distance of closest approach to rational surface is DEL (read from TJ JSON file)
  int            NFIX;    // Number of fixups (read from TJ JSON file)
  int            NDIAG;   // Number of radial grid-points for diagnostics (read from TJ JSON file)
  double         NULC;    // Use zero pressure jump conditions when |nu_L| < NULC (read from TJ JSON file)
  int            ITERMAX; // Maximum number of iterations used to determine quantities at rational surface (read from TJ JSON file)

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

  double*            cmuw;       // cosh(mu) on wall
  double*            eetaw;      // eta on wall
  double*            cetaw;      // cos(eta) on wall
  double*            setaw;      // sin(eta) on wall
  double*            R2grgzw;    // R^2 nabla r . nabla z   on wall
  double*            R2grgew;    // R^2 nabla r . nabla eta on wall

  gsl_spline*        Rrzwspline; // Interpolated R2grgz function on wall
  gsl_spline*        Rrewspline; // Interpolated R2grge function on wall
  gsl_spline*        Rwspline;   // Interpolated R function on wall
  gsl_spline*        Zwspline;   // Interpolated Z function on wall

  gsl_interp_accel*  Rrzwacc;    // Accelerator for interpolated R2grgz function on wall
  gsl_interp_accel*  Rrewacc;    // Accelerator for interpolated R2grge function on wall
  gsl_interp_accel*  Rwacc;      // Accelerator for interpolated R function on wall
  gsl_interp_accel*  Zwacc;      // Accelerator for interpolated Z function on wall

  // -----------
  // Island data
  // -----------
  int                Nh;         // Number of harmonics in island calculation (read from ISLAND JSON file)
  int                NX;         // Number radial grid-points in island calculation (read from ISLAND JSON file)
  double*            T0inf;      // Asymptotic value of X - deltaT[0] at given rational surface (read from Outputs/Island/Island.nc)
  double*            XX;         // Island solution radial grid (read from Outputs/Island/Island.nc)
  Array<double,3>    deltaTh;    // Harmonics of perturbed electron temperature in vicinity of island at given rational surface (read from Outputs/Island/Island.nc)

  gsl_spline**       dThspline;  // Interpolated island harmonic functions
  gsl_interp_accel** dThacc;     // Accelerator for interpolated island harmonic functions
  
  // --------------------------
  // Plasma boundary parameters
  // --------------------------
  double* tbound; // theta values on plasma boundary
  double* Rbound; // R values on plasma boundary
  double* Zbound; // Z values on plasma boundary
  double* dRdthe; // dR/dtheta values on plasma boundary
  double* dZdthe; // dZ/dtheta values on plasma boundary
  double  igrr2b; // <|nabla r|^(-2)> on plasma boundary

  // ---------------
  // Wall parameters
  // ---------------
  double* twall;  // theta values on wall
  double* Rwall;  // R values on wall
  double* Zwall;  // Z values on wall
  double* dRdthw; // dR/dtheta values on wall
  double* dZdthw; // dZ/dtheta values on wall
  double  H1w;    // H1 at wall
  double  H1pw;   // H1p at wall
  double  igrr2w; // <|nabla r|^(-2)> on plasma boundary
  double  alphaw; // Wall effective radius parameter
  double  gammaw; // Normalized thin-wall resistive wall mode growth-rate

  // -----------------------------------
  // Resonant magnetic perturbation data
  // -----------------------------------
  int     ncoil; // Number of toroidal strands that make up RMP coils 
  double* Rcoil; // R coodinates of strands 
  double* Zcoil; // Z coodinates of strands 
  double* Icoil; // Toroidal currents flowing in strands 

  // --------------------
  // Vacuum solution data
  // --------------------
  double                   sa;     // Edge magnetic shear

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
  double                   bw;     // Relative wall radius (read from Equilibrium JSON file)

  Array<complex<double>,2> Pwal;   // Wall solution matrix
  Array<complex<double>,2> Qwal;   // Wall solution matrix
  Array<complex<double>,2> Rwal;   // Wall solution matrix
  Array<complex<double>,2> Swal;   // Wall solution matrix

  Array<complex<double>,2> iImat;  // Rwal * iImat = Swal
  Array<complex<double>,2> iIher;  // Hermitian component of iImat
  Array<complex<double>,2> iIant;  // Anti-Hermitian component of iImat

  Array<complex<double>,2> PImat;  // PImat = Pvac * iImat - Qvac
  Array<complex<double>,2> RImat;  // RImat = Rvac * iImat - Svac
  Array<complex<double>,2> IRmat;  // IRmat = Rvac * iImat 

  Array<complex<double>,2> RPImat; // RPImat * RImat = PImat
  Array<complex<double>,2> RPIdag; // Hermitian conjugate of RPImat
  Array<complex<double>,2> Gmat;   // Perfect-wall vacuum response matrix: Gmat = (1/2) (RPImat + RPIdag)

  Array<complex<double>,2> iRPImat; // iRPImat * PImat = RImat
  Array<complex<double>,2> iRPIdag; // Hermitian conjugate of iRPImat
  Array<complex<double>,2> iGmat;   // Perfect-wall vacuum response matrix: iGmat = (1/2) (iRPImat + iRPIdag)

  Array<complex<double>,2> Rbamat;  // Rwal Rvac^-1
   
  // ---------------------
  // Rational surface data
  // ---------------------
  int     nres;    // Number of rational magnetic flux-surfaces
  int*    mres;    // Poloidal mode numbers at rational surfaces
  double* qres;    // Safety-factors at rational surfaces
  double* rres;    // Minor radii of rational surfaces
  double* Pres;    // PsiN values at rational surfaces
  double* qerr;    // Residual in determination of rational surface
  double* sres;    // Magnetic shears at rational surfaces
  double* gres;    // g values at rational surfaces
  double* hres;    // h values at rational surfaces
  double* DIres;   // DI values at rational surfaces
  double* DRres;   // DR values at rational surfaces
  double* nuLres;  // Ideal Mercier indices of large solution at rational surfaces
  double* nuSres;  // Ideal Mercier indices of small solution at rational surfaces
  int*    Jres;    // Index of resonant poloidal harmonic at rational surfaces
  double* Flarge;  // TJ/STRIDE scaling factors for large solution
  double* Fsmall;  // TJ/STRIDE scaling factors for small solution
  double* neres;   // Electron number densities at rational surfaces
  double* nepres;  // Radial gradients of electron number densities at rational surfaces
  double* Teres;   // Electron temperatures at rational surfaces
  double* Tepres;  // Radial gradients of electron temperatures at rational surfaces
  double* Ls;      // Magnetic shear-lengths at rational surfaces
  double* LT;      // Pressure gradient scale-lengths at rational surfaces
  double* Lc;      // Average magnetic field-line curvature scale-lengths at rational surfaces
  double* betah;   // Normalized electron pressure at rational surfaces
  double* alphab;  // Bootstrap parameters at rational surfaces
  double* alphac;  // Curvature parameters at rational surfaces
 
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
  double*                  Pgrid; // PsiN values for grid-points for diagnostics
  Array<double,2>          Ttest; // Torque test for solution vectors versus radius
  Array<double,2>          Pnorm; // Norms of Psi components of solution vectors versus radius
  Array<double,2>          Znorm; // Norms of Z components of solution vectors versus radius
  double*                  hode;  // Step-length versus radius
  double*                  eode;  // Truncation error versus radius
  Array<complex<double>,3> YYY;   // Solution vectors versus radius
  Array<complex<double>,2> Pi;    // Reconnected fluxes at rational surfaces associated with solution vectors
  Array<complex<double>,1> dPi;   // Current sheets at rational surfaces associated with solution vectors

  // -------------------------------------
  // Tearing mode dispersion relation data
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
  Array<complex<double>,3> psiu;  // Scaled psi components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> zu;    // z components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> chiu;  // chi components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> xiu;   // xi components of unreconnected tearing eigenfunctions
  Array<double,2>          Tf;    // Torques associated with fully reconnected eigenfunctions
  Array<double,2>          Tu;    // Torques associated with unreconnected eigenfunctions
  Array<double,3>          Tfull; // Torques associated with pairs of fully reconnected eigenfunctions
  Array<double,3>          Tunrc; // Torques associated with pairs of unreconnected eigenfunctions

  // ----------------------------
  // Temperature and density data
  // ----------------------------
  double*                  Psik;  // Reconnected fluxes at rational surfaces
  complex<double>*         PsTp;  // Reconnected fluxes for temperature calculation outside rational surfaces
  complex<double>*         PsTm;  // Reconnected fluxes for temperature calculation inside rational surfaces
  double*                  dTp;   // Temperature adjustments outside rational surfaces
  double*                  dTm;   // Temperature adjustments inside rational surfaces
  complex<double>*         Psnp;  // Reconnected fluxes for electron number density calculation outside rational surfaces
  complex<double>*         Psnm;  // Reconnected fluxes for electron number density calculation inside rational surfaces
  double*                  dnp;   // Electron number density adjustments outside rational surfaces
  double*                  dnm;   // Electron number density  adjustments inside rational surfaces
  double*                  delta; // Island asymmetry parameters at rational surfaces
  double*                  width; // Island widths at rational surfaces
  
  Array<complex<double>,3> neu;   // n_e components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> Teu;   // T_e components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> dneu;  // delta n_e components of unreconnected tearing eigenfunctions
  Array<complex<double>,3> dTeu;  // delta T_e components of unreconnected tearing eigenfunctions

  // -----------------------------------
  // Resonant magnetic perturbation data
  // -----------------------------------
  complex<double>*         Psix;    // RMP perturbation at plasma boundary
  complex<double>*         Xi;      // RMP response vector
  complex<double>*         Upsilon; // RMP response vector
  complex<double>*         Lambda;  // RMP response vector
  complex<double>*         Chi;     // RMP drive at rational surfaces
  Array<complex<double>,2> Psirmp;  // Psi component of ideal RMP response eigenfunction
  Array<complex<double>,2> Zrmp;    // Z component of ideal RMP response eigenfunction
  complex<double>*         Psixs;   // Psix on plasma boundary
  complex<double>*         Psirmps; // Psirmp on plasma boundary

  // ---------------------------------------
  // Visualization of tearing eigenfunctions
  // ---------------------------------------
  int                      Nf;     // Number of radial grid-points on visualization grid (from Equilibrium.nc)
  int                      Nw;     // Number of poloidal grid-points on visualization grid (from Equilibrium.nc)
  double*                  rf;     // Radial grid-points on visualization grid (from Equilibrium.nc)
  Array<double,2>          RR;     // R coordinates of visualization grid-points (from Equilibrium.nc)
  Array<double,2>          ZZ;     // Z coordinates of visualization grid-points (from Equilibrium.nc)
  Array<double,2>          dRdr;   // dR/dr values at visualization grid-points (from Equilibrium.nc)
  Array<double,2>          dRdt;   // dR/dtheta values at visualization grid-points (from Equilibrium.nc)
  Array<double,2>          dZdr;   // dZ/dr values at visualization grid-points (from Equilibrium.nc)
  Array<double,2>          dZdt;   // dZ/dtheta values at visualization grid-points (from Equilibrium.nc)
  Array<double,2>          rvals;  // r values of visualization grid-points (from Equilibrium.nc)
  Array<double,2>          thvals; // theta values of visualization grid-points (from Equilibrium.nc)
  Array<double,2>          Bmod;   // |B| on visualization grid
  Array<double,2>          Btor;   // |B_toroidal| on visualization grid
  Array<double,2>          Bpol;   // |B_poloidal| on visualization grid
  Array<complex<double>,3> Psiuf;  // Psi components of Fourier-transformed unreconnected tearing eigenfunctions 
  Array<complex<double>,3> Zuf;    // Z components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> psiuf;  // Scaled psi components of Fourier-transformed unreconnected tearing eigenfunctions 
  Array<complex<double>,3> zuf;    // z components of Fourier-transformed unreconnected tearing eigenfunctions 
  Array<complex<double>,3> chiuf;  // chi components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> xiuf;   // xi^r components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> neuf;   // n_e components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> Teuf;   // T_e components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> dneuf;  // delta n_e components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> dTeuf;  // delta T_e components of Fourier-transformed unreconnected tearing eigenfunctions
  Array<complex<double>,3> Psiuv;  // Psi components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,3> Zuv;    // Z components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,3> zuv;    // z components of unreconnected tearing eigenfunctions on visualization grid
  Array<complex<double>,3> chiuv;  // Chi components of unreconnected tearing eigenfunctions on visualization grid
  Array<double,3>          bRc;    // Cosine component of b^R on visualization grid
  Array<double,3>          bRs;    // Sin component of b^R on visualization grid
  Array<double,3>          bZc;    // Cosine component of b^Z on visualization grid
  Array<double,3>          bZs;    // Sin component of b^Z on visualization grid
  Array<double,3>          bPc;    // Cosine component of R b^phi on visualization grid
  Array<double,3>          bPs;    // Sin component of R b^^phi on visualization grid
  Array<double,3>          xic;    // Cosine component of xi^r on visualization grid
  Array<double,3>          xis;    // Sin component of xi^r on visualization grid

  // ----------------------------------------
  // Visualization of temperature and density 
  // ----------------------------------------
  int                      NPHI;   // Number of toroidal grid-points on extended vizualization grid (from TJ JSON file)
  double*                  PP;     // Toroidal grid-points on extended visualization grid
  Array<double,4>          nec;    // n_e on extended visualization grid
  Array<double,4>          Tec;    // T_e on extended visualization grid
  Array<double,4>          dnec;   // delta n_e on extended visualization grid
  Array<double,4>          dTec;   // delta T_e on extended visualization grid

  // -------------------------------------
  // Visualization of ideal eigenfunctions
  // -------------------------------------
  Array<complex<double>,3> Xief;   // Xi components of Fourier transformed no-wall ideal eigenfunctions on visualization grid
  Array<complex<double>,3> pXief;  // Xi components of Fourier transformed perfect-wall ideal eigenfunctions on visualization grid
  Array<complex<double>,3> Xiev;   // Xi components of no-wall ideal eigenfunctions on visualization grid
  Array<complex<double>,3> pXiev;  // Xi components of perfect-wall ideal eigenfunctions on visualization grid

  // -----------------------------
  // Synthetic ECE diagnostic data 
  // -----------------------------
  double  tilt;           // Tilt angle of central chord (degrees) (read from Equilibrium JSON file)
  double* req;            // r values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* weq;            // omega values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* teq;            // theta values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* Req;            // R values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* Zeq;            // Z values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* BReq;           // Equilibrium B_parallel values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* neeq;           // Equilibrium n_e values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* Teeq;           // Equilibrium T_e values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* dRdreq;         // dR/dr values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* dRdteq;         // (dR/dt)/r values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* dZdreq;         // dZ/dr values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* dZdteq;         // (dZ/dt)/r values on central chord (read from Outputs/Equilibrium/Equilibrium.nc)
  double* Leq;            // Length along central chord, L
  double* Lres;           // L values of rational surfaces on central chord
  double* Rres;           // R values of rational surfaces on central chord 
  double* Ores;           // R values of island O-points on central chord
  double* Xres;           // R values of island X-points on central chord

  Array<double,3> bReqc;  // Perturbed B_parallel values on central chord
  Array<double,3> neeqc;  // Total n_e values on central chord
  Array<double,3> Teeqc;  // Total T_e values on central chord
  Array<double,3> dneeqc; // Perturbed n_e values on central chord
  Array<double,3> dTeeqc; // Perturbed T_e values on central chord

  double* DeltaO;  // Inward radial shift of 1st harmonic O-mode convolution function on central chord
  double* sigmaO;  // Standard deviation of 1st harmonic O-mode convolution function on central chord
  double* tauO;    // Optical depth of 1st harmonic O-mode convolution function on central chord
  double* DeltaX;  // Inward radial shift of 2nd harmonic X-mode convolution function on central chord
  double* sigmaX;  // Standard deviation of 2nd harmonic X-mode convolution function on central chord
  double* tauX;    // Optical depth of 2nd harmonic X-mode convolution function on central chord

  gsl_spline**       Teeqcspline;  // Interpolated total T_e values on central chord
  gsl_spline**       dTeeqcspline; // Interpolated perturbed T_e values on central chord
  gsl_interp_accel** Teeqcacc;     // Accelerator for interpolated total T_e values on central chord
  gsl_interp_accel** dTeeqcacc;    // Accelerator for interpolated perturbed T_e values on central chord
  Array<double,3>    TeeqdO;       // Total T_e values on central chord modified by relativistic 1st harmonic O-mode ece broadening
  Array<double,3>    dTeeqdO;      // Perturbed T_e values on central chord modified by relativistic 1st harmonic O-mode ece broadening
  Array<double,3>    TeeqdX;       // Total T_e values on central chord modified by relativistic 2nd harmonic X-mode ece broadening
  Array<double,3>    dTeeqdX;      // Perturbed T_e values on central chord modified by relativistic 2nd harmonic X-mode ece broadening

  // --------------------------------------------------------
  // Visualization of resonant magnetic perturbation response
  // --------------------------------------------------------
  Array<complex<double>,2> Psirf;  // Psi components of Fourier-transformed ideal RMP response eigenfunction
  Array<complex<double>,2> Zrf;    // Z components of Fourier-transformed ideal RMP response eigenfunction
  Array<complex<double>,2> Psirv;  // Psi components of ideal RMP response eigenfunction
  Array<complex<double>,2> Zrv;    // Z components of ideal RMP response eigenfunction

  // ----------------------------
  // No-wall ideal stability data
  // ----------------------------
  Array<complex<double>,3> Psii;    // Psi components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> Zi;      // Z components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> Xii;     // Xi components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> xii;     // xi components of ideal solutions launched from magnetic axis
  Array<complex<double>,3> Chii;    // Chi components of ideal solutions launched from magnetic axis
  Array<complex<double>,2> Ji;      // Poloidal harmonics of current on plasma boundary associated with ideal solutions launched from magnetic axis
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
  Array<complex<double>,3> Psie;    // Psi components of ideal eigenfunctions
  Array<complex<double>,3> Ze;      // Z components of ideal eigenfunctions
  Array<complex<double>,3> Xie;     // Xi components of ideal eigenfunctions
  Array<complex<double>,3> xie;     // xi components of ideal eigenfunctions
  Array<complex<double>,2> Je;      // Poloidal harmomics of current on plasma boundary associated with ideal eigenfunctions
  double*                  deltaW;  // Total perturbed ideal potential energy 
  double*                  deltaWp; // Plasma contribution to perturbed ideal potential energy
  double*                  deltaWv; // Vacuum contribution to perturbed ideal potential energy
  Array<complex<double>,2> Psiy;    // Psi values on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> Jy;      // Current on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> Xiy;     // Xi values on plasma boundary associated with ideal eigenfunctions
  complex<double>*         gammax;  // Expansion of Psi_x at boundary in ideal eigenfunctions
  complex<double>*         gamma;   // Expansion of Psi_rmp at boundary in ideal eigenfunctions
  complex<double>*         ya;      // psi values of zeroth eigenfunction at plasma boundary
  int                      jzero;   // Index of first solution with positive vacuum energy

  // ------------------------------
  // Perfect-wall ideal energy data
  // ------------------------------
  Array<complex<double>,2> pWmat;    // Plasma ideal energy matrix
  Array<complex<double>,2> pWher;    // Hermitian component of Wmat
  Array<complex<double>,2> pWant;    // Anti-Hermitian component of Wmat
  double*                  pWval;    // Eigenvalues of symmeterized W-matrix
  Array<complex<double>,2> pVmat;    // Vacuum ideal energy matrix
  Array<complex<double>,2> pVher;    // Hermitian component of Vmat
  Array<complex<double>,2> pVant;    // Anti-Hermitian component of Vmat
  double*                  pVval;    // Eigenvalues of symmeterized V-matrix
  Array<complex<double>,2> pUmat;    // Total ideal energy matrix
  Array<complex<double>,2> pUher;    // Hermitian component of Umat
  Array<complex<double>,2> pUant;    // Anti-Hermitian component of Umat
  double*                  pUval;    // Eigenvalues of symmeterized U-matrix
  Array<complex<double>,2> pUvec;    // Eigenvectors of symmeterized U-matrix
  Array<complex<double>,2> pUres;    // Residuals of Uvec orthonormality matrix
  Array<complex<double>,3> pPsie;    // y components of ideal eigenfunctions
  Array<complex<double>,3> pZe;      // Z components of ideal eigenfunctions
  Array<complex<double>,3> pXie;     // Xi components of ideal eigenfunctions
  Array<complex<double>,3> pxie;     // xi components of ideal eigenfunctions
  Array<complex<double>,2> pJe;      // Poloidal harmomics of current on plasma boundary associated with ideal eigenfunctions
  double*                  pdeltaW;  // Total perturbed ideal potential energy 
  double*                  pdeltaWp; // Plasma contribution to perturbed ideal potential energy
  double*                  pdeltaWv; // Vacuum contribution to perturbed ideal potential energy
  Array<complex<double>,2> pPsiy;    // y values on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> pJy;      // Current on plasma boundary associated with ideal eigenfunctions
  Array<complex<double>,2> pXiy;     // Z values on plasma boundary associated with ideal eigenfunctions
  complex<double>*         yb;       // psi values of zeroth eigenfunction at wall
  int                      pjzero;   // Index of first solution with positive vacuum energy
 
  // -------------------
  // Resonant layer data
  // -------------------
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
  double* Wdres;  // Critical island width for temperature flattening

  // ----
  // Misc
  // ----
  int    rhs_chooser, iomega;
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

  // ......................
  // In ReadEquilibrium.cpp
  // ......................

  // Read equilibrium data from Outputs/Equilibrium/Equilibrium.nc
  void ReadEquilibrium ();
  // Calculate metric data at plasma boundary
  void CalculateMetricBoundary ();
  // Calculate metric data at wall
  void CalculateMetricWall ();
  // Read island data from Outputs/Island/Island.nc
  void ReadIsland ();
    
  // .............
  // In Vacuum.cpp
  // .............

  // Calculate vacuum boundary matrices
  void GetVacuumBoundary ();
  // Calculate vacuum wall matrices
  void GetVacuumWall ();
  // Evaluate right-hand sides of vacuum odes
  void CashKarp45Rhs1 (double r, complex<double>* Y, complex<double>* dYdr) override;
  // Evaluate eta for RMP coil calculation
  double Geteta (double R, double Z, double Rp, double Zp);
  // Evaluate G for RMP coil calculation
  double GetG (double R, double Z, double Rp, double Zp);
  
  // ...............
  // In Rational.cpp
  // ...............
  // Target function for finding rational surfaces
  double RootFindF (double r) override;
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
  // Get value of kmp 
  double Getkmp (double r, int m);
  // Get value of km 
  double Getkm (double r, int m);
  // Get value of Lmm
  double GetLmm (double r, int m);
  
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

  // ..........
  // In RMP.cpp
  // ..........

  // Calculate resonant magnetic perturbation data
  void CalculateResonantMagneticPerturbation ();

  // ................
  // In Visualize.cpp
  // ................
  
  // Calculate unreconnected eigenfunction visualization data
  void VisualizeEigenfunctions ();
  // Calculate resonant magnetic response visualization data
  void VisualizeRMP ();
  // Evaluate right-hand sides of ece odes
  void CashKarp45Rhs (double R, double* Y, double* dYdr) override;

  // ............
  // In Ideal.cpp
  // ............

  // Calculate no-wall ideal stability
  void CalculateNoWallIdealStability ();
  // Calculate perfect-wall ideal stability
  void CalculatePerfectWallIdealStability ();

  // ..........
  // In RWM.cpp
  // ..........

  // Calculate resistive wall mode stability
  void CalculateRWMStability ();
   
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
  // Return value of P1
  double GetP1 (double r);
  // Return value of P2
  double GetP2 (double r);
  // Return value of P3
  double GetP3 (double r);
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
  // Return value of DI
  double GetDI (double r);
  // Return value of DR
  double GetDR (double r);
  // Return value of Flarge
  double GetFlarge (double r, int m);
  // Return value of Fsmall
  double GetFsmall (double r, int m);
  // Return value of Ls
  double GetLs (double r);
  // Return value of LT
  double GetLT (double r);
  // Return value of Lc
  double GetLc (double r);
  // Return value of betae
  double Getbetah (double r);
  // Return value of alphab
  double Getalphab (double r);
  // Return value of alphac
  double Getalphac (double r);

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
  // In Netcdf.cpp
  // .............

  // Read equilibrium data from netcdf file
  void ReadEquilibriumNetcdf ();
  // Read island data from netcdf file
  void ReadIslandNetcdf (int k);
  // Write stability data to netcdf file
  void WriteNetcdf ();
};
