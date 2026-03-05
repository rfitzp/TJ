// ########################################################################################################
// Program to calculate stability boundaries and satuarted island widths for modified Wesson-profile plasma
// ########################################################################################################

#pragma once

#include "Utility.h"

// ############
// Class header
// ############
class Wesson : private Utility
{
 private:
  
  // ------------------
  // Physics parameters
  // ------------------
  double mpol; // Poloidal mode number (read from JSON file)
  double ntor; // Toroidal mode number  (read from JSON file)
  double rw;   // Radius of wall (relative to minor radius of plasma) (read from JSON file)
  double q0;   // Safety-factor at magnetic axis (read from JSON file)
  double qa;   // Safety factor at edge of plasma (read from JSON file)
  double abar; // Radius of current channel (read from JSON file)
  double qafc; // End of qa scan is qafc * qs (read from JSON file)

  int    flag; // Control flag  file)
               // flag = 0 - stability limit scan
               // flag = 1 - 2/1 Delta scan at q0=0.8
  
  // ----------------------
  // Integration parameters
  // ----------------------
  double eps;   // Closest distance to magnetic axis
  double del;   // Closest distance to rational surface
  double h0;    // Initial step-length
  double acc;   // Integration accuracy
  double hmin;  // Minimum step-length

  // ---------------------------------------
  // One-dimensional root finding parameters
  // ---------------------------------------
  int    Nint;    // Number of search intervals
  double Eta;     // Minimum magnitude of F at root F(x) = 0
  int    Maxiter; // Maximum number of iterations

  // ----
  // Misc
  // ----
  int count, flg;
  
  // ----------------------
  // Public class functions
  // ----------------------
 public:
  
  Wesson ();              // Constructor
  virtual ~Wesson () {};  // Destructor

  void Solve  ();         // Solve problem

  // -----------------------
  // Private class functions
  // -----------------------
 private:

  // Calculate tearing stability index
  double GetDelta ();
  // Find resonant surface radius
  double Findrs (double qs);
  // Get value of parameter nu
  double Getnu ();
  // Get current profile
  double GetJ (double r);
  // Get first derivative of current profile
  double GetJp (double r);
  // Get second derivaitive of current profile
  double GetJpp (double r);
  // Get safety-factor profile
  double Getq (double r);
  // Get magnetic shear profile
  double Gets (double r);
  // Get value of parameter alpha
  double Getalpha (double r);
  // Get value of parameter beta
  double Getbeta (double r);
  // Get value of normalized self-inductance
  double Getli ();
  // Get value of normalized self-inductance
  double Getli1 ();

  // Evaluate right-hand sides of differential equations
  void Rhs (double x, Array<double,1>& y, Array<double,1>& dydx);
  // Adaptive step integration routine
  void RK4Adaptive (double& x, Array<double,1>& y, double& h, double& t_err, 
		    double acc, double S, int& rept, int maxrept, 
		    double h_min, double h_max, int flag, int diag, FILE* file);
  // Fixed step integration routine
  void RK4Fixed (double& x, Array<double,1>& y, double h);

  // Evaluate target function for one-dimensional root finding
  double RootFindF (double x);
  // One-dimensional root finding routine
  double RootFind (double x1, double x2);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double x1, double x2, double F1, double F2, double& x);

  FILE* OpenFile (char* filename);
};

