// Layer.h

// #####################################################################

// Class to solve three-field resonant layer equations in tokamak plasma

// Inputs:
//  Inputs/TJ.json          - TJ JSON file
//  Inputs/Layer.json       - Layer JSON file

// Outputs:
//  Outputs/Layer/Layer.nc

// Plotting scripts:
//  Plots/Layer/*.py

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

// Documentation: ../Documentation/Layer.pdf

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
class Layer : private Utility
{
private:

  // .............
  // TJ parameters
  // .............
  double* input;      // TJ input parameters
  int     nres;       // Number of rational surfaces
  double* r_res;      // Radii of rational surfaces (read from Outputs/TJ/TJ.nc)
  int*    m_res;      // Poloidal mode numbers of rational surfaces (read from Outputs/TJ/TJ.nc)
  double* S13_res;    // Cube roots of Lundquist numbers at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* tau_res;    // Normalized resistive kink timescales at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* QE_res;     // Normalized ExB frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Qe_res;     // Normalized electron diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Qi_res;     // Normalized ion diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* iotae_res;  // Ratio of minus electron diamagnetic frequencies to total diamagnetic frequencies at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* D_res;      // Normalized ion sound radii at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Pphi_res;   // Normalized momentum diffusivities at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Pperp_res;  // Normalized energy diffusivities at rational surfaces (read from Outputs/TJ/TJ.nc
  double* Delta_res;  // Delta_primes at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Deltac_res; // Critical Delta_primes at rational surfaces (read from Outputs/TJ/TJ.nc)
  double* Chi_res;    // Absolute values of Chi (read from Outputs/TJ/TJ.nc)

  // .........................
  // Layer equation parameters
  // .........................
  double          Qe;     // Normalized electron diamagnetic frequency 
  double          Qi;     // Normalized ion diamagnetic frequency
  double          QE;     // Normalized ExB frequency
  double          D;      // Normalized ion sound radius 
  double          Pphi;   // Normalized momentum diffusivity 
  double          Pperp;  // Normalized energy diffusivity 
  double          iotae;  // Ratio of minus electron diamagnetic frequency to total diamagnetic frequency

  double          g_r;    // Real part of normalized growth-rate in MHD frame
  double          g_i;    // Imaginary part of normalized growth-rate in MHD frame
  double          Delta;  // Normalized effective tearing stability index
  complex<double> Deltas; // Normalized complex layer response index

  // ....................................
  // Layer calculation control parameters
  // ....................................
  int    RMP;    // Flag to enable resonant magnetic perturbation calculation (read from TJ JSON file)
  int    MARG;   // Flag to enable calculation of ion branch marginal stability points (read from JSON file)
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
  double* f_e;     // Electron-branch relative frequency (+1/0/-1 if mode corotates with electron/ExB/ion fluid)
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

  // ....
  // Misc
  // ....
  int lowD;
  complex<double> Im; // Square-root of -1: Im = complex<double> (0., 1.)
  
public:

  // Constructor
  Layer ();

  // Solve problem
  void Solve (int verbose);

private:

  // Read TJ data from netcdf file
  void ReadNetcdf ();
  // Find marginal stability points associated with ith rational surface
  void FindMarginal (int i, int verbose);
  // Find complex growth-rate of electron branch-mode associated with ith rational surface
  void GetElectronBranchGrowth (int i, int verbose);
  // Find shielding factor and locking torque curves associated with ith rational surface
  void GetTorque (int i);
  // Write Layer data to netcdf file
  void WriteNetcdf ();
  // Deallocate memory
  void CleanUp ();

  // Solve layer equations 
  void SolveLayerEquations ();
  
  // Evaluate right-hand sides of differential equations
  void CashKarp45Rhs (double x, complex<double>*  y, complex<double>*  dydx) override;

  // Target functions for two-dimensional root finding via Newton-Raphson method
  void NewtonFunction (double x1, double x2, double& F1, double& F2) override;

  // Target function for one-dimensional root finding
  double RootFindF (double x) override;
};
