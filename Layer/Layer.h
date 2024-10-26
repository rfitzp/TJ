// Layer.h

// #####################################################################

// Class to solve three-field resonant layer equations in tokamak plasma

// Inputs:
//  Inputs/Layer.json - JSON file

// Outputs:
//  Outputs/Layer/Layer.nc

// Plotting scripts:
//  Plots/Layer/*.py

// Class uses following external libraries:
//  Blitz++ library        (https://github.com/blitzpp/blitz)
//  nclohmann JSON library (https://github.com/nlohmann/json)
//  netcdf-c++ library     (https://github.com/Unidata/netcdf-cxx4)

// Author:
// Richard Fitzpatrick,
// Institute of Fusion Studies,
// Department of Physics
// University of Texas at Austin
// rfitzp@utexas.edu

// Source: https://github.com/rfitzp/TJ

// Documentation: ../Documentation/Layer.pdf

// #####################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <blitz/array.h>
#include <nlohmann/json.hpp>
#include <netcdf>

using namespace blitz;
using           json = nlohmann::json;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

// ############
// Class header
// ############
class Layer
{
private:

  // .............
  // TJ parameters
  // .............
  double* input;      // TJ input parameters
  int     nres;       // Number of rational surfaces
  double* r_res;      // Radii of rational surfaces (read from Outputs/TJ/TJ.nc)
  int*    m_res;      // Poloidal mode numbers of rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Delta_res;  // Delta_primes at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Deltac_res; // Critical Delta_primes at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Chi_res;    // Absolute values of Chi (read from Outputs/TJ/TJ.nc)
  double* S13_res;    // Cube roots of Lundquist numbers at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* tau_res;    // Normalized resistive kink timescales at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* QE_res;     // Normalized ExB frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Qe_res;     // Normalized electron diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Qi_res;     // Normalized ion diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* iotae_res;  // Ratio of electron diamagnetic frequencies to total diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* D_res;      // Normalized ion sound radii at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Pphi_res;   // Normalized momentum diffusivities at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Pperp_res;  // Normalized energy diffusivities at rational surfaces (read from Outputs/TJ/TJ.nc

  // .........................
  // Layer equation parameters
  // .........................
  double          Qe;     // Normalized electron diamagnetic frequency 
  double          Qi;     // Normalized ion diamagnetic frequency 
  double          D;      // Normalized ion sound radius 
  double          Pphi;   // Normalized momentum diffusivity 
  double          Pperp;  // Normalized energy diffusivity 
  double          iotae;  // Ratio of electron diamagnetic frequency to total diamagnetic frequency

  double          g_r;    // Real part of normalized growth-rate in MHD frame
  double          g_i;    // Imaginary part of normalized growth-rate in MHD frame
  double          Delta;  // Normalized effective tearing stability index
  complex<double> Deltas; // Normalized complex layer response index

  // ....................................
  // Layer calculation control parameters
  // ....................................
  int    LAYER;  // Switch to enable layer calculation (read from JSON file)
  double pstart; // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double pend;   // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double P3max;  // Value of Pmax[3] above which switch to low-D layer equations made (read from JSON file)
  int    Nscan;  // Number of points in marginal stability and frequency scans (read from JSON file)

  // .........................
  // Marginal stability points
  // .........................
  Array<int,1>    np_marg;  // Number of marginal stability points at rational surfaces 
  Array<double,2> gr_marg;  // Real part of normalized growth-rates in MHD frame at marginal stability points at rational surfaces 
  Array<double,2> gi_marg;  // Imaginary part of normalized growth-rates in MHD frame at marginal stability points at rational surfaces 
  Array<double,2> Dr_marg;  // Real part of normalized tearing stability indicesat marginal stability points at rational surfaces 
  Array<double,2> Di_marg;  // Imaginary part of normalized tearing stability indices at marginal stability points at rationaL surfaces
  
  // .................................
  // Growth-rates and real frequencies
  // .................................
  double* gamma_e; // Electron-branch growth-rates at rational surfaces (kHz)
  double* omega_e; // Electron-branch real frequencies at rational surfaces  (kHz)
  double* res_e;   // Electron-branch residuas of zero search at rational surfaces 
  int*    lowD_e;  // Electron-branch low D flags at rational surfaces

  // ..........................................
  // Shielding factor and locking torque curves
  // ..........................................
  Array<double,2> omega_r; // Real frequencies of RMPs
  Array<double,2> Deltar;  // Real part of Delta_layer
  Array<double,2> Deltai;  // Imaginary part of Delta_layer
  Array<double,2> Xi_res;  // Shielding factor curves
  Array<double,2> T_res;   // Locking torque curves

  // ...............................
  // Adaptive integration parameters
  // ...............................
  double acc;     // Integration accuracy (read from JSON file)
  double h0;      // Initial integration step-length (read from JSON file)
  double hmin;    // Minimum integration step-length (read from JSON file)
  double hmax;    // Maximum integration step-length (read from JSON file)
  int    maxrept; // Maximum number of step recalculations
  int    flag;    // Error calculation flag

  // ...........................
  // Newton iteration parameters
  // ...........................
  double eps;     // Step-size for calculation of Jacobian (read from JSON file)
  double smax;    // Maximum step-size (read from JSON file)
  double smin;    // Minimum step-size (read from JSON file)
  double Eta;     // Minimum magnitude of f at root f(x) = 0 (read from JSON file)
  int    Maxiter; // Maximum number of iterations (read from JSON file)

  // .........................
  // Cash-Karp RK45 parameters
  // .........................
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ....
  // Misc
  // ....
  int count, lowD;
  complex<double> Im; // Square-root of -1: Im = complex<double> (0., 1.)
  
public:

  // Constructor
  Layer ();

  // Solve problem
  void Solve ();

private:

  // Read TJ data from netcdf file
  void ReadNetcdf ();
  // Find marginal stability points associated with ith rational surface
  void FindMarginal (int i);
  // Find complex growth-rate of electron branch-mode associated with ith rational surface
  void GetElectronBranchGrowth (int i);
  // Find shielding factor and locking torque curves associated with ith rational surface
  void GetTorque (int i);
  // Write Layer data to netcdf file
  void WriteNetcdf ();
  // Deallocate memory
  void CleanUp ();

  // Solve layer equations 
  void SolveLayerEquations ();
  
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, complex<double>*  y, complex<double>*  dydx);
  // Advance set of coupled first-order o.d.e.s by single step using adaptive step-length
  //  Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Adaptive (int neqns, double& x, complex<double>* y, double& h, 
			   double& t_err, double acc, double S, double T, int& rept,
			   int maxrept, double h_min, double h_max, int flag);
  // Advance set of coupled first-order o.d.e.s by single step using fixed step-length
  //  Cash-Karp fourth-order/fifth-order Runge-Kutta scheme
  void CashKarp45Fixed (int neqns, double& x, complex<double>* y, complex<double>* err, double h);

  // Get Jacobian matrix
  void GetJacobian (double& J11, double& J12, double& J21, double& J22);
  // Find root of Deltas = Delta via Newton iteration
  void GetRoot (double& gr, double& gi, double& F);
  // Target function for zero finding
  double Feval (double x);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double x1, double x2, double F1, double F2, double& x);

  // Strip comments from a string
  string stripComments (const string& input);
  // Read JSON file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (const char* filename);
};
