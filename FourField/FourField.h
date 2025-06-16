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

#include <blitz/array.h>
#include <netcdf>
#include "Utility.h"

using namespace blitz;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

// ############
// Class header
// ############
class FourField : private Utility
{
private:

  // ..............................
  // Four field equation parameters
  // ..............................
  double          g_r;     // Real part of normalized growth-rate in MHD frame (read from JSON file)
  double          g_i;     // Imaginary part of normalized growth-rate in MHD frame (read from JSON file)
  double          Qe;      // Normalized electron diamagnetic frequency (read from JSON file)
  double          Qi;      // Normalized ion diamagnetic frequency (read from JSON file)
  double          D;       // Normalized ion sound radius (read from JSON file)
  double          Pphi;    // Normalized momentum diffusivity (read from JSON file)
  double          Pperp;   // Normalized energy diffusivity (read from JSON file)
  double          cbeta;   // Normalized plasma pressure (read from JSON file)
  
  double          iotae;   // Ratio of minus electron diamagnetic frequency to total diamagnetic frequency: -Qe /(Qi - Qe)

  complex<double> Deltas3; // Normalized three-field complex layer response index
  complex<double> Deltas4; // Normalized four-field complex layer response index
  complex<double> Deltas5; // Normalized four-field complex layer response index

  // ..........................................
  // Four field  calculation control parameters
  // ..........................................
  double pstart;  // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double pend;    // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double P3max;   // Value of Pmax[3] above which switch to low-D three-field layer equations made (read from JSON file)

  // ....
  // Misc
  // ....
  int lowD, rhs_chooser;
  complex<double> Im; // Square-root of -1: Im = complex<double> (0., 1.)
  
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
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, complex<double>*  y, complex<double>*  dydx) override;
};
