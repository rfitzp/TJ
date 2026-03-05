// NonLinear.h

// ################################################################################
// Class to calculate nonlinear coupling of three tearing modes in a toroidal pinch
// ################################################################################

// Parallel current profile:
//
//  sigma(r) = (2*epsa/q0) (1 - r^alphas)^nus
//
// Pressure profile:
//
//  P(r) = (beta0/2) (1 - r^alphap)^nup
// 

#ifndef NONLINEAR
#define NONLINEAR

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

// ############
// Class header
// ############
class NonLinear
{
 private:

  // ----------------------
  // Equilibrium parameters
  // ----------------------
  int    machine; // 1 - RFP; 2 - Tokamak
  double epsa;    // Inverse aspect-ratio of plasma boundary
  double epsb;    // Inverse aspect-ratio of perfectly conducting wall
  double q0;      // Central safety-factor
  double beta0;   // Central plasma beta
  double alphas;  // Parallel current profile parameter
  double nus;     // Parallel current profile parameter
  double alphap;  // Pressure profile parameter
  double nup;     // Pressure profile parameter

  double Theta;   // Pinch parameter (zero pressure)
  double TTheta;  // Pinch parameter (continuous pressure)
  double Frev;    // Reversal parameter (zero pressure)
  double FFrev;   // Reversal parameter (continuous pressure)
  
  // ---------------
  // Mode parameters
  // ---------------
  double mpol;   // Poloidal mode number
  double ntor;   // Toroidal mode number
  double mpol1;  // Poloidal mode number of first mode
  double ntor1;  // Toroidal mode number of first mode
  double mpol2;  // Poloidal mode number of second mode
  double ntor2;  // Toroidal mode number of second mode
  double mpol3;  // Poloidal mode number of third mode
  double ntor3;  // Toroidal mode number of third mode

  double rs1;    // Rational surface radius of first mode
  double rs2;    // Rational surface radius of second mode
  double rs3;    // Rational surface radius of third mode

  double W1;     // Island width at rational surface of first mode
  double W2;     // Island width at rational surface of second mode
  double W3;     // Island width at rational surface of third mode

  double Delta1; // Tearing stabilty index of first mode (no pressure jumps)
  double Delta2; // Tearing stabilty index of second mode (no pressure jumps)
  double Delta3; // Tearing stabilty index of third mode (no pressure jumps)

  double* Psi1;  // Tearing eigenfuction of first mode (no pressure jumps)
  double* Psip1; // Tearing eigenfuction derivative of first mode (no pressure jumps)
  double* Psi2;  // Tearing eigenfuction of second mode (no pressure jumps)
  double* Psip2; // Tearing eigenfuction derivative of second mode (no pressure jumps)
  double* Psi3;  // Tearing eigenfuction of third mode (no pressure jumps)
  double* Psip3; // Tearing eigenfuction derivative of third mode (no pressure jumps)

  double* JJ1;   // J values at control surfaces for first mode
  double* JJ2;   // J values at control surfaces for second mode
  double* JJ3;   // J values at control surfaces for third mode

  double DDelta1; // Tearing stabilty index of first mode (with pressure jumps)
  double DDelta2; // Tearing stabilty index of second mode (with pressure jumps)
  double DDelta3; // Tearing stabilty index of third mode (with pressure jumps)

  double* PPsi1;  // Tearing eigenfuction of first mode (with pressure jumps)
  double* PPsip1; // Tearing eigenfuction derivative of first mode (with pressure jumps)
  double* PPsi2;  // Tearing eigenfuction of second mode (with pressure jumps)
  double* PPsip2; // Tearing eigenfuction derivative of second mode (with pressure jumps)
  double* PPsi3;  // Tearing eigenfuction of third mode (with pressure jumps)
  double* PPsip3; // Tearing eigenfuction derivative of third mode (with pressure jumps)

  double* Xi1;    // Regualarized displacement of first mode (with pressure jumps)
  double* Xi2;    // Regualarized displacement of second mode (with pressure jumps)
  double* Xi3;    // Regualarized displacement of third mode (with pressure jumps)

  // ----------------------
  // Calculation parameters
  // ----------------------
  int    Ngrid;  // Number of radial grid points
  int    Ncntr;  // Number of control surfaces
  double eps;    // Closest approach to magnetic axis
  double delta;  // Closest approach to rational surface
  int    TorCur; // Toroidal curvature flag:
                 // 0 - no toroidal curvature correction
                 // 1 - correct jump conditions for magnetic field
                 // 2 - correct jump conditions for eigenfunction
   
  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;       // Integration accuracy
  double h0;        // Initial step-length
  double hmin;      // Minimum step-length
  double hmax;      // Maximum step-length
  int    maxrept;   // Maximum number of step recalculations
  int    flag;      // Integration error calculation flag

  // ----------------
  // Calculation data
  // ----------------
  double* rr;       // Radial grid
  double* ssigma;   // Parallel current profile (continuous pressure)
  double* PP;       // Pressure profile (continuous pressure)
  double* PQ;       // Modified pressure profile (continuous pressure)
  double* BBphi;    // Toroidal magnetic field (continuous pressure)
  double* BBtheta;  // Poloidal magnetic field (continuous pressure)
  double* qq;       // Safety-factor profile (continuous pressure)
  double* qp;       // Derivative of safety-factor profile (continuous pressure)

  gsl_interp_accel* APQ; // Interpolation accelerator for modified pressure profile
  gsl_spline*       SPQ; // Interpolator for modified pressure profile

  double* rrc;      // Radii of control surfaces
  double* PPc;      // Pressures between control surfaces
  double* Btheta2j; // Jump in square of poloidal field across control surfaces
  double* Bthetam;  // Poloidal magnetic fields to left of control surfaces
  double* Bthetap;  // Poloidal magnetic fields to right of control surfaces

  double* sssigma;  // Parallel current profile (stepped pressure)
  double* PPP;      // Pressure profile (stepped pressure)
  double* BBBphi;   // Toroidal magnetic field (stepped pressure)
  double* BBBtheta; // Poloidal magnetic field (stepped pressure)
  double* qqq;      // Safety-factor profile (stepped pressure)

  // ------------------
  // Overlap parameters
  // ------------------
  double*           Tau;    // Overlap function (no pressure jumps)
  double*           Taup;   // Overlap function (with pressure jumps)
  gsl_interp_accel* Atau;   // Interpolation accelerator for overlap function (no pressure jumps)
  gsl_interp_accel* Ataup;  // Interpolation accelerator for overlap function (with pressure jumps)
  gsl_spline*       Stau;   // Interpolator for overlap function (no pressure jumps)
  gsl_spline*       Staup;  // Interpolator for overlap function (with pressure jumps)
  double*           Over;   // Overlap integral (no pressure jumps)
  double*           Overp;  // Overlap integral (with pressure jumps)

  // ------------------
  // RK4/RK5 parameters
  // ------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int count, rhs_chooser;
  
public:

  // Constructor
  NonLinear (int _machine);
  // Destructor
  ~NonLinear ();

  // Set beta0
  void Setbeta0 (double _beta0);
  // Set mode numbers
  void SetMode (int _n, int _k);
  // Set W
  void SetW (double _W);
  // Solve problem
  void Solve (int verbose, double& tau, double& taup, double& _Delta1, double& _Delta2, double& _Delta3);
  // Evaluate stability index
  void GetDelta (int _mpol, int _ntor, double& Delta, double& Psirv);

private:

   // Calculate continuous zero pressure equilibrium
  void CalcZero ();
  // Calculate continuous pressure equilibrium
  void CalcEquilibrium ();
  // Calculate stepped pressure equilibrium
  void CalcStepped ();
  // Calculate overlap function
  void GetTau ();
  // Calculate overlap integral
  void CalcOverlap ();

  // Calculate tearing eigenfunction (with no pressure jumps)
  void GetPsi (double _mpol, double _ntor, double rs, double& Delta, double* Psi, double* Psip, int verbose);
  // Calculate tearing eigenfunction (with pressure jumps)
  void GetPsiStepped (double _mpol, double _ntor, double rs, double W, double* JJ, double& DDelta, double* PPsi, double* PPsip, double* Xi, int verbose);
  // Calculate tearing stability index (with pressure jumps)
  void GetDelta (double _mpol, double _ntor, double rs, double* JJ, double& Delta, double& Psirv);
  // Find rational surface radius
  double Findrs (double qs);

  // Get F
  double GetF (double r, double _mpol, double _ntor);
  // Get regularized 1/F
  double GetinvF (double r, double rs, double W, double _mpol, double _ntor);
  // Get G
  double GetG (double r, double _mpol, double _ntor);
  // Get H
  double GetH (double r, double _mpol, double _ntor);
  
  // Get J values at control surfaces
  void GetJ (double _mpol, double _ntor, double rs, double W, double* JJ);
  // Get K
  double GetK (double r, double rs, double W);
  // Get 1 /(m - n q)
  double Getimnq (double r, double rs, double W);
  // Get f
  double Getf (double r);
  // Get g
  double Getg (double r);
  // Get lambda
  double Getlambda (double r);
  // Get F
  double GetF (double r);
  // Get Fp
  double GetFp (double r);
  // Get G
  double GetG (double r);
  // Get H
  double GetH (double r);
  // Get toroidal magnetic field
  double GetBphi (double r);
  // Get poloidal magnetic field
  double GetBtheta (double r);
  // Get safety-factor
  double Getq (double r);
  // Get derivative of safety-factor
  double Getqp (double r);
  
  // Get equilibrium parallel current
  double GetSigma (double r);
  // Get equilibrium parallel current gradient
  double GetSigmap (double r);
  // Get equilibrium pressure
  double GetP (double r);
  // Get equilibrium pressure gradient
  double GetPp (double r);
  // Get equilibrium pressure
  double GetPQ (double r);
  // Get equilibrium pressure gradient
  double GetPQp (double r);

  // Evaluate right-hand sides of differential equations
  void Rhs (double x, double* y, double* dydx);

  // Adaptive step length RK4/RK5 integration routine
  void RK4RK5Adaptive (int neqns, double& x, double* y, double& h, 
		       double& t_err, double acc, double S, double T, int& rept,
		       int maxrept, double h_min, double h_max, int flag, 
		       int diag, FILE* file);
  // Fixed step length RK4/RK5 integration routine
  void RK4RK5Fixed (int neqns, double& x, double* y, double* err, double h);

  // Open new file for writing
  FILE* OpenFilew (char* filename);
};

#endif //NONLINEAR
