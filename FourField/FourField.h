// FourField.h

// #####################################################################

// Class to solve four-field resonant layer equations in tokamak plasma

// Inputs:
//  Inputs/FourField.json       - FourField JSON file

// Outputs:
//  Outputs/FourField/FourField.nc

// Plotting scripts:
//  Plots/FourField/*.py

// Class uses following external libraries:
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)
//  GNU scientific library (https://www.gnu.org/software/gsl)

// Author:
// Richard Fitzpatrick,
// Institute of Fusion Studies,
// Department of Physics
// University of Texas at Austin
// rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/FourField.pdf

// #####################################################################

#pragma once

#include <gsl/gsl_sf_gamma.h>

#include "Utility.h"

// ############
// Class header
// ############
class FourField : private Utility
{
private:

  // ..............................
  // Four-field equation parameters
  // ..............................
  double          g_r;     // Real part of normalized growth-rate in MHD frame (read from JSON file)
  double          g_i;     // Imaginary part of normalized growth-rate in MHD frame (read from JSON file)
  double          Qe;      // Normalized electron diamagnetic frequency (read from JSON file)
  double          Qi;      // Normalized ion diamagnetic frequency (read from JSON file)
  double          D;       // Normalized ion sound radius (read from JSON file)
  double          Pphi;    // Normalized momentum diffusivity (read from JSON file)
  double          Pperp;   // Normalized energy diffusivity (read from JSON file)
  double          cbeta;   // Normalized plasma pressure (read from JSON file)

  int             JK;      // Flag for J.-K. Park calculation in which Pperp = cbeta^2 (read from JSON file)
  
  double          iotae;   // Ratio of minus electron diamagnetic frequency to total diamagnetic frequency: -Qe /(Qi - Qe)

  complex<double> Deltas3; // Normalized three-field complex layer response index
  complex<double> Deltas4; // Normalized four-field complex layer response index
  complex<double> Deltas5; // Normalized four-field complex layer response index

  // .........................................
  // Four-field calculation control parameters
  // .........................................
  double pstart;  // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double pend;    // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double P3max;   // Value of Pmax[3] above which switch to low-D three-field layer equations made (read from JSON file)

  // ....................................
  // Four-field frequency scan parameters
  // ....................................
  int    Scan;    // Flag for performing frequency scan (read from JSON file)
  double Fscan;   // Scan from g_i = - Fmax * (Qe - Qi) to g_i = Fmax * (Qe - Qi) (read from JSON file)
  int    Nscan;   // Number of points in scan (read from JSON file)

  double* givals;  // g_i values in frequency scan
  double* D3rvals; // Re(Delta3) values in frequency scan
  double* D3ivals; // Im(Delta3) values in frequency scan
  double* D4rvals; // Re(Delta4) values in frequency scan
  double* D4ivals; // Im(Delta4) values in frequency scan

  double* J1rvals; // Re(DeltaJ1) values in frequency scan
  double* J1ivals; // Im(DeltaJ1) values in frequency scan
  double* J2rvals; // Re(DeltaJ2) values in frequency scan
  double* J2ivals; // Im(DeltaJ2) values in frequency scan

  // ....
  // Misc
  // ....
  int lowD, rhs_chooser;
  
public:

  // Constructor
  FourField ();

  // Solve problem
  void Solve ();

private:

  // Write layer data to netcdf file
  void WriteNetcdf ();
  // Deallocate memory
  void CleanUp ();

  // Solve three-field layer equations 
  void SolveThreeFieldLayerEquations ();
  // Solve four-field layer equations 
  void SolveFourFieldLayerEquations ();
  // Calculate Jace Waybright's approximate Deltas
  void GetJace (complex<double>& _DeltaJ1, complex<double>& _DeltaJ2);
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, complex<double>*  y, complex<double>*  dydx) override;
};
