// ECE.h

// #####################################################################################

// Class to calculate spatial convolution function for perpendicular 1st harmonic O-mode
// and 2nd harmonic X-mode electron cyclotron emission

// Inputs:
//  Inputs/ECE.json - JSON file

// Outputs:
//  Outputs/ECE/ECE.nc

// Plotting scripts:
//  Plots/ECE/*.py

// Class uses following external libraries:
//  nclohmann JSON library             (https://github.com/nlohmann/json)
//  GNU scientific library             (https://www.gnu.org/software/gsl)
//  Steven Johsnson's Faddeeva pacakge (http://ab-initio.mit.edu/faddeeva)
//  netcdf-c++ library                 (https://github.com/Unidata/netcdf-cxx4)

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

// ############
// Class header
// ############
class ECE : private Utility
{
 private:

  // ------------------
  // Physics parameters
  // ------------------
  double Te;  // Electron temperature at cyclotron resonance (eV) (read from JSON file)
  double ne;  // Electron number density at cyclotron resonance (m^-3) (read from JSON file)
  double B0;  // Toroidal magnetic field-strength on magnetic axis (T) (read from JSON file)
  double R0;  // Major radius of magnetic axis (m) (read from JSON file)
  double Rw;  // Major radius of cyclotron resonance (m) (read from JSON file)

  // ------------------
  // Derived parameters
  // ------------------
  double wc0;  // Cyclotron frequency at magnetic axis (s^-1)
  double vt;   // Electron thermal speed at cyclotron resonance (m/s)
  double vtc2; // (v_t/c)^2 at cyclotron resonance
  double tau0; // Optical depth parameter
  double wp;   // Electron plasma frequency at cyclotron resonance (m/s)
  double w1;   // Ratio of emission frequency to plasma frequency
  double tw;   // Ratio of electron thermal energy to rest mass energy

  // ----------------------
  // Calculation parameters
  // ----------------------
  double  zmax;   // Maximum value of z in z-array (read from JSON file)
  int     znum;   // Number of points in z-array (read from JSON file)
  double* zz;     // z-array
  double* F72r;   // Real part of F72 evaluated on z-array
  double* F72i;   // Imaginary part of F72 evaluated on z-array

  double  wcmin;  // Minumum value of wc/w1 in wc-array (read from JSON file)
  double  wcmax;  // Maximum value of wc/w1 in wc-array (read from JSON file)
  int     wcnum;  // Number of points in wc-array (read from JSON file)
  double* wwc;    // wc-array
  double* alphaO; // 1st harmonic O-mode absorption coefficient evaluated on wc-array
  double* alphaX; // 2nd harmonic X-mode absorption coefficient evaluated on wc-array
  double* tauO;   // 1st harmonic O-mode optical depth evaluated on wc-array
  double* tauX;   // 2nd harmonic X-mode optical depth evaluated on wc-array
  double* HO;     // 1st harmonic O-mode spectral convolution function evaluated on wc-array
  double* HX;     // 2nd harmonic X-mode spectral convolution function evaluated on wc-array
  double* RR;     // Major radii evaluated on wc-array
  double* FO;     // 1st harmonic O-mode spatial convolution function evaluated on wc-array
  double* FX;     // 2nd harmonic X-mode spatial convolution function evaluated on wc-array

  double sigmaO;  // Standard deviation of 1st harmonic O-mode spatial convolution function
  double sigmaX;  // Standard deviation of 2nd harmonic X-mode spatial convolution function
  double DeltaO;  // Inward radial shift of 1st harmonic O-mode spatial convolution function
  double DeltaX;  // Inward radial shift of 2nd harmonic X-mode spatial convolution function
  double tauOst;  // Saturated optical depth of 1st harmonic O-mode
  double tauXst;  // Saturated optical depth of 2nd harmonic X-mode

  double* FOfit;  // Truncated Gaussian fit to 1st harmonic O-mode spatial convolution function evaluated on wc-array
  double* FXfit;  // Truncated Gaussian fit to 2nd harmonic X-mode spatial convolution function evaluated on wc-array

  // ----
  // Misc
  // ----
  int rhs_chooser;
   
 public:

  // Constructor
  ECE ();
  // Destructor
  ~ECE ();

  // Solve problem
  void Solve ();
  // Calculate standard deviation and radial shift of 1st harmonic O-mode spatial convolution function
  void GetOmodeFit (double _Te, double _ne, double _B0, double _R0, double _Rw, double& sigma, double& Delta, double& taust);
  // Calculate standard deviation and radial shift of 2nd harmonic X-mode spatial convolution function
  void GetXmodeFit (double _Te, double _ne, double _B0, double _R0, double _Rw, double& sigma, double& Delta, double& taust);

 private:
  
  // Set derived calculation parameters
  void SetDerivedParameters (double _Te, double _ne, double _B0, double _R0, double _Rw);

  // Evaluate ECE absorption function F_7/2
  complex<double> Get_F72 (double z);
  // Get resonance function
  double Get_z (double wc);
  // Get refractive index of 1st harmonic O-mode
  double Get_NperpO (double wc);
  // Get refractive index of 2nd harmonic X-mode
  double Get_NperpX (double wc);
  // Get absorption coefficient of 1st harmonic O-mode
  double Get_alpha1O (double wc);
  // Get absorption coefficient of 2nd harmonic X-mode
  double Get_alpha2X (double wc);

  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, double*  y, double*  dydx) override;

  // Write netcdf file
  void Write_netcdf ();
};
