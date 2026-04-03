// BetaLimit.h

// Program to search for beta limit as function of central safeety-factor using TJ code

#pragma once

#include "Utility.h"
#include "TJ.h"

// ############
// Class header
// ############
class BetaLimit : private Utility
{
 private:

  double qc;        // Central safety-factor
  double rs;        // Rational surface radius
  double ss;        // Resonant magnetic shear
  double dW;        // delta W
  
  double qc_start;  // Lower limit of central safety-factor scan (read from JSON file)
  double qc_end;    // Upper limit of central safety-factor scan (read from JSON file)
  int    N_qc;      // Number of points in central safety-factor scan (read from JSON file)

  double pc_start;  // Lower limit of central pressure search interval (read from JSON file)
  double pc_end  ;  // Upper limit of central pressure search interval (read from JSON file)

 public:
  
  // Constructor
  BetaLimit ();
  // Destructor
  ~BetaLimit ();

  // Solve problem
  void Solve ();

 private:

  // Target function for one-dimensional root finding
  double RootFindF (double x) override;
};
