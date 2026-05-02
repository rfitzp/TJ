// BetaLimit.h

// Program to search for beta limit as function of specified parameter using TJ code

#pragma once

#include "Utility.h"
#include "TJ.h"

// ############
// Class header
// ############
class BetaLimit : private Utility
{
 private:

  double val;        // Value of scanned parameter
  double rs;         // Rational surface radius
  double ss;         // Resonant magnetic shear
  double dW;         // delta W

  int    CNTRL;      // Sets which parameter is scanned (read from JSON file)
                     // 0 - qc scan; 1 - bw scan; 2 - Ea scan
  int    UP;         // Sets whether beta limit goes up (1) or down (0) with increasing parameter (read from JSON file)
  double val_start;  // Lower limit of parameter scan (read from JSON file)
  double val_end;    // Upper limit of parameter scan (read from JSON file)
  int    N_val;      // Number of points in parameter scan (read from JSON file)
 
  double pc_start;   // Lower limit of central pressure search interval (read from JSON file)
  double pc_end  ;   // Upper limit of central pressure search interval (read from JSON file)

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
