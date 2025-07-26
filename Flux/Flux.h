// Flux.h

// ####################################################################################
// Class to read EFIT file and generate TJ coordinate system in order to check accuracy
// of EFIT data

// .................
// Calculation grid:
// .................

// Radial grid in PsiN = 1. - Psi /Psi_axis (assuming Psi = 0 on boundary) is 
//
// PsiN_j = PSIVAC * [1. - (1. - s)^PACK]^2 for j = 0, NPSI-1
//
//  where s = (j+1) /NPSI
//
// Poloidal grid in PEST poloidal angle theta is
//
// theta_k = 2.*M_PI * t for k = 0, NTHETA-1
//
//  where t = k /(NTHETA-1).
//
// theta = 0 on inboard midplane.
// theta > 0 above midplane.
//
// All lengths normalized to R_0 (major radius of magnetic axis).
// All magnetic field-strengths normalized to B_0 (on-axis vacuum toroidal magnetic field).

// ...................
// Inputs and outputs:
// ...................
// Calculation control parameters in JSON file Inputs/Flux.json

// EFIT file in Outputs/WriteEFIT/EFIT.txt
// Data written to Outputs/Flux

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// #####################################################################################

#pragma once

#include "Utility.h"

// Pointers to right-hand side function for adaptive integration
extern "C" int pRhs1 (double, const double[], double[], void*);
extern "C" int pRhs2 (double, const double[], double[], void*);
extern "C" int pRhs3 (double, const double[], double[], void*);

// gFile reading function
extern "C" void gFileRead ();

// ############
// Class header
// ############
class Flux : private Utility
{
 private:

  // .............................................
  // Control parameters read from Inputs/Flux.json
  // .............................................
  
  int    NPSI;    // Number of points in PsiN grid
  double PACK;    // Packing index for PsiN grid
  int    NTHETA;  // Number of points in theta grid

  double H0;      // Initial integration step-length for equilibrium flux-surface integrals 
  double ACC;     // Integration accuracy for equilibrium flux-surface integrals

  // ..................
  // Stage 1 parameters
  // ..................
  
  double  R0;     // Scale major radius (m)
  double  B0;     // Scale toroidal magnetic field strength (T)
  double  RLEFT;  // Bounding box coordinate 
  double  ZLOW;   // Bounding box coordinate
  double  RRIGHT; // Bounding box coordinate
  double  ZHIGH;  // Bounding box coordinate

  double  Raxis;  // Magnetic axis coordinate
  double  Zaxis;  // Magnetic axis coordinate
  int     NRPTS;  // Number of R points
  double* RPTS;   // R array
  int     NZPTS;  // Number of Z points
  double* ZPTS;   // Z array
  int     NBPTS;  // Number of boundary points
  double* RBPTS;  // R on plasma boundary
  double* ZBPTS;  // Z on plasma boundary
  int     NLPTS;  // Number of limiter points
  double* RLPTS;  // R on limiter boundary
  double* ZLPTS;  // Z on limiter boundary
  
  double* PSIN;   // PsiN array 
  double* G;      // g(Psi)
  double* Pr;     // p(Psi) 
  double* GGp;    // g dg/dPsi 
  double* Prp;    // dp/dPsi
  double* Q;      // q(Psi)

  Array<double,2> PSIARRAY; // Psi (R, Z) array (unnormalized)
  double*         PSIAXIS;  // Interpolated Psi array at Z-level of magnetic axis

  // .......................................
  // Flux coordinate construction parameters
  // .......................................
  
  double  Psic;   // Psi on magnetic axis
  double  Rbound; // R coordinate of plasma boundary on inboard mid-plane
  int     ia;     // R grid index of inboard plasma boundary at Z = Z_axis
  int     ic;     // R grid index of magnetic axis
  int     jc;     // Z grid index of magnetic axis
  int     L;      // Number of points in Psi(R,Zaxis) array
  double* s;      // Array of s = sqrt[1 - Psi(R,Zaxis)] values
  double* Rs;     // Array of R(s) values, where R is major radius on inboard midplane
  double  qa;     // Safety-factor at plasma boundary
  double  ra;     // Radial coordinate of plasma boundary
  double  qgp;    // Parameter passed to integration right-hand sides

  // .........................
  // Stage2 profile parameters
  // ..........................
  
  double* P;    // Psi array
  double* RP;   // R(Psi) on inboard mid-plane
  double* rP;   // r(Psi)
  double* GP;   // g(Psi)
  double* QGP;  // q(Psi)/g(Psi)
  double* QP;   // q(Psi)
  double* FP;   // f(Psi)
  double* PP;   // P(Psi)
  double* GPP;  // dg/dPsi
  double* PPP;  // dP/dPsi
  double* GPX;  // g dg/dPsi (check)
  double* PPX;  // dP/dPsi (check)
  double* S;    // sqrt(1 - Psi)
  double* QX;   // q(Psi) from EFIT file
  double* PsiN; // PsiN array

  double*     th;   // theta array
  gsl_matrix* RRst; // R versus theta on general surfaces
  gsl_matrix* ZZst; // Z versus theta on general surfaces
  gsl_matrix* RRr;  // (dR/dr)_theta
  gsl_matrix* RRth; // (dR/dtheta)_r /r
  gsl_matrix* ZZr;  // (dZ/dr)_theta
  gsl_matrix* ZZth; // (dZ/dtheta)_r /r
  gsl_matrix* Jac;  // Jacobian
  gsl_matrix* Jax;  // Analytic Jacobian

  // ...........................
  // Plasma boundary data for TJ
  // ...........................
  
  double* Rbndry;    // R values on boundary
  double* Zbndry;    // Z values on boundary
  double* dRdtheta;  // dRdtheta values on boundary
  double* dZdtheta;  // dZdtheta values on boundary

 public:

  // ...........
  // In Flux.cpp
  // ...........
  
  // Constructor
  Flux ();
  // Destructor
  ~Flux ();

  // Solve problem
  void Solve ();

  // ................
  // In Integrate.cpp
  // ................
  
  // Evaluate right-hand sides of q/g equation
  int Rhs1 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of r equation
  int Rhs2 (double r, const double y[], double dydr[], void*);
  // Evaluate right-hand sides of theta equation
  int Rhs3 (double r, const double y[], double dydr[], void*);
 
private:

  // ...........
  // In Flux.cpp
  // ...........

  // Set global parameters
  void SetParameters ();
  // Input gFile data and output Stage1 data
  void Stage1 ();

  // .............
  // In Stage2.cpp
  // .............

  // Input Stage1 data and output Stage2 data
  void Stage2 ();
  // Write Stage2 NETCDF file
  void WriteStage2Netcdfc ();

  // .............
  // In Plasma.cpp
  // .............
  
  // Read data for Stage2 calculations
  void Stage2ReadData ();
  // Calculate Stage2 q profile
  void Stage2CalcQ ();
  // Calculate Stage2 straight angle coordinate system
  void Stage2CalcStraightAngle ();
  // Calculate boundary data for TJ
  void Stage2CalcBoundary ();

  // ................
  // In Integrate.cpp
  // ................
 
  // Calculate q(P)/g(P) profile
  void CalcQGP ();
  // Calculate r(P) profile
  void CalcrP ();
  // Calculate straight angle data on general magnetic flux-surfaces
  void CalcStraightAngleGeneral ();

  // ..................
  // In Interpolate.cpp
  // ..................
  
  // 1D interpolation function with nonuniform grid
  double Interpolate         (int I, double* X, double* Y, double x,                                 int order);
  double InterpolateCubic    (       double* X, double* Y, double x, int i0, int i1, int i2,         int order);
  double InterpolateQuartic  (       double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order);

  double Interpolate1        (int I, double* X, gsl_matrix* Y, int i, double x,                                 int order);
  double InterpolateCubic1   (       double* X, gsl_matrix* Y, int i, double x, int i0, int i1, int i2,         int order);
  double InterpolateQuartic1 (       double* X, gsl_matrix* Y, int i, double x, int i0, int i1, int i2, int i3, int order);

  // Interpolate on uniform 2D grid
  double Interpolate               (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX,
				    int order, int cubic);
  double InterpolateCubicCubic     (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX,
				    int i0, int i1, int i2, int j0, int j1, int j2,                 int order);
  double InterpolateQuarticCubic   (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX,
				    int i0, int i1, int i2, int i3, int j0, int j1, int j2,         int order);
  double InterpolateCubicQuartic   (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX,
				    int i0, int i1, int i2, int j0, int j1, int j2, int j3,         int order);
  double InterpolateQuarticQuartic (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX,
				    int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order);

  // Evaluate Psi (R, Z)
  double GetPsi  (double r, double z);
  // Evaluate dPsi/dR (R, Z)
  double GetPsiR (double r, double z);
  // Evaluate dPsi/dZ (R, Z)
  double GetPsiZ (double r, double z);
};

