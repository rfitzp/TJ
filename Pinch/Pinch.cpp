// Pinch.cpp

#include "Pinch.h"

#define NINPUT 9

// ###########
// Constructor
// ###########
Pinch::Pinch ()
{
  // ---------------------------------------------
  // Ensure that directory ../Outputs/Pinch exists
  // ---------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Pinch"))
    {
      exit (1);
    }

  // ......................................
  // Read control parameters from JSON file
  // ......................................
  string JSONFilename = "../Inputs/Pinch.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  epsa   = JSONData["epsa"]  .get<double> ();
  q0     = JSONData["q0"]    .get<double> ();
  beta0  = JSONData["beta0"] .get<double> ();
  alphas = JSONData["alphas"].get<double> ();
  nus    = JSONData["nus"]   .get<double> ();
  alphap = JSONData["alphap"].get<double> ();
  nup    = JSONData["nup"]   .get<double> ();

  Ngrid = JSONData["Ngrid"]  .get<int>    ();
  eps   = JSONData["eps"]    .get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (epsa <= 0.)
    {
      printf ("Pinch:: Error - epsa must be positive\n");
      exit (1);
    }
  if (q0 <= 0.)
    {
      printf ("Pinch:: Error - q0 must be positive\n");
      exit (1);
    }
  if (beta0 < 0.)
    {
      printf ("Pinch:: Error - beta0 cannot be negative\n");
      exit (1);
    }
  if (alphas < 0.)
    {
      printf ("Pinch:: Error - alphas cannot be negative\n");
      exit (1);
    }
  if (nus < 0.)
    {
      printf ("Pinch:: Error - nus cannot be negative\n");
      exit (1);
    }
  if (alphap < 0.)
    {
      printf ("Pinch:: Error - alphap cannot be negative\n");
      exit (1);
    }
  if (nup < 0.)
    {
      printf ("Pinch:: Error - nup cannot be negative\n");
      exit (1);
    }
  if (Ngrid < 1)
    {
      printf ("Pinch:: Error - Ngrid must be larger than unity\n");
      exit (1);
    }
  if (eps < 0. || eps > 1.)
    {
      printf ("Pinch:: Error - eps must lie in the interval 0 to 1\n");
      exit (1);
    }
  
  printf ("\n");
  printf ("Class PINCH::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("epsa  = %10.3e q0  = %10.3e beta0 = %10.3e alphas = %10.3e nus = %10.3e alphap = %10.3e nup = %10.3e\n",
	  epsa, q0, beta0, alphas, nus, alphap, nup);
  printf ("Ngrid =  %4d      eps = %10.3e\n",
	  Ngrid, eps);
  
  // ---------------
  // Allocate memory
  // ---------------
  rr       = new double[Ngrid+1];
  ssigma   = new double[Ngrid+1];
  PP       = new double[Ngrid+1];
  BBphi    = new double[Ngrid+1];
  BBtheta  = new double[Ngrid+1];
  qq       = new double[Ngrid+1];

  qqc      = new double[Ngrid+1];
  PPc      = new double[Ngrid+1];
  BBphic   = new double[Ngrid+1];
  BBthetac = new double[Ngrid+1];

  PPp      = new double[Ngrid+1];
  PPpc     = new double[Ngrid+1];
  PPpm     = new double[Ngrid+1];

  // ------------------
  // Set up radial grid
  // ------------------
  double dr = 1. /double (Ngrid);
  for (int i = 0; i <= Ngrid; i++)
    rr[i] = double (i) * dr;
}

// ###########
// Destructor
// ###########
Pinch::~Pinch ()
{
  delete[] rr;  delete[] ssigma; delete[] PP;     delete[] BBphi;    delete[] BBtheta; delete[] qq;
  delete[] qqc; delete[] PPc;    delete[] BBphic; delete[] BBthetac;
  delete[] PPp; delete[] PPpc;   delete[] PPpm;
}

// #########################
// Function to solve problem
// #########################
void Pinch::Solve ()
{
  // Calculate toroidal pinch equilibrium
  CalcEquilibrium ();
 
  // Output data to netcdf file
  WriteNetcdf ();
}

// ######################################
// Function to output data to netcdf file
// ######################################
void Pinch::WriteNetcdf ()
{
   printf ("Writing data to netcdf file Outputs/Pinch/Pinch.nc:\n");

   double Input[NINPUT];
   
   Input[0] = epsa;
   Input[1] = q0;
   Input[2] = beta0;
   Input[3] = alphas;
   Input[4] = nus;
   Input[5] = alphap;
   Input[6] = nup;
   Input[7] = double (Ngrid);
   Input[8] = eps;
      
   try
     {
       NcFile dataFile ("../Outputs/Pinch/Pinch.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim i_d = dataFile.addDim ("Ni", NINPUT);
       NcDim r_d = dataFile.addDim ("Nr", Ngrid+1);
         
       NcVar i_x     = dataFile.addVar ("InputParameters", ncDouble, i_d);
       i_x.putVar (Input);
       NcVar r_x     = dataFile.addVar ("r",               ncDouble, r_d);
       r_x.putVar (rr);
       NcVar sigma_x = dataFile.addVar ("sigma",           ncDouble, r_d);
       sigma_x.putVar (ssigma);
       NcVar P_x     = dataFile.addVar ("P",               ncDouble, r_d);
       P_x.putVar (PP);
       NcVar Bt_x    = dataFile.addVar ("B_phi",           ncDouble, r_d);
       Bt_x.putVar (BBphi);
       NcVar Bp_x    = dataFile.addVar ("B_theta",         ncDouble, r_d);
       Bp_x.putVar (BBtheta);
       NcVar q_x     = dataFile.addVar ("q",               ncDouble, r_d);
       q_x.putVar (qq);
       NcVar Pc_x    = dataFile.addVar ("P_c",             ncDouble, r_d);
       Pc_x.putVar (PPc);
       NcVar Btc_x   = dataFile.addVar ("B_phi_c",         ncDouble, r_d);
       Btc_x.putVar (BBphic);
       NcVar Bpc_x   = dataFile.addVar ("B_theta_c",       ncDouble, r_d);
       Bpc_x.putVar (BBthetac);
       NcVar qc_x    = dataFile.addVar ("q_c",             ncDouble, r_d);
       qc_x.putVar (qqc);
       NcVar pp_x    = dataFile.addVar ("dPdr",            ncDouble, r_d);
       pp_x.putVar (PPp);
       NcVar ppc_x   = dataFile.addVar ("dPdr_crit",       ncDouble, r_d);
       ppc_x.putVar (PPpc);
       NcVar ppm_x   = dataFile.addVar ("dPdr_marg",       ncDouble, r_d);
       ppm_x.putVar (PPpm);
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/Pinch/Pinch.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }
}

// #####################
// Function to set beta0
// #####################
void Pinch::Setbeta0 (double _beta0)
{
  beta0 = _beta0;
}

// #####################
// Function to set q0
// #####################
void Pinch::Setq0 (double _q0)
{
  q0 = _q0;
}

// ###################
// Function to set nus
// ###################
void Pinch::Setnus (double _nus)
{
  nus = _nus;
}

// #################################
// Function to calculate equilibrium
// #################################
void Pinch::CalcEquilibrium ()
{
  // ................................................
  // Calculate parallel current and pressure profiles
  // ................................................
  for (int i = 0; i <= Ngrid; i++)
    {
      ssigma[i] = GetSigma (rr[i]);
      PP    [i] = GetP     (rr[i]);
    }

  // ...................................................
  // Calculate magnetic field and safety-factor profiles
  // ...................................................
  rhs_chooser = 0;
  
  BBphi[0]   = 1.;
  BBtheta[0] = 0.;
  qq[0]      = q0;
 
  double  r, h, t_err;
  int     rept;
  double* y   = new double[3];
  double* err = new double[3];

  r     = eps;
  h     = h0;
  count = 0;

  y[0] = 1.;
  y[1] = (epsa /q0) * r;
  y[2] = 0.;

  for (int i = 1; i <= Ngrid; i++)
    {
      do
	{
	  CashKarp45Adaptive (3, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (3, r, y, err, rr[i] - r);

      BBphi[i]   = y[0];
      BBtheta[i] = y[1];
      qq[i]      = epsa * rr[i] * BBphi[i] /BBtheta[i];
     }

  Theta = y[1] /2./y[2];
  Frev  = y[0] /2./y[2];
  
  printf ("\nF = %10.4e Theta = %10.4e\n", Frev, Theta);
  
  delete[] y; delete[] err;

  // .............................................................................
  // Calculate Mercier-stable magnetic field, safety-factor, and pressure profiles
  // .............................................................................
  rhs_chooser = 1;
  
  BBphic[0]   = 1.;
  BBthetac[0] = 0.;
  qqc[0]      = q0;
  PPc[0]      = beta0/2.;

  PPp[0]      = 0.;
  PPpc[0]     = 0.;
  PPpm[0]     = 0.;
 
  y   = new double[4];
  err = new double[4];

  r     = eps;
  h     = h0;
  count = 0;

  y[0] = 1.;
  y[1] = (epsa /q0) * r;
  y[2] = 0.;
  y[3] = beta0/2.;

  for (int i = 1; i <= Ngrid; i++)
    {
      do
	{
	  CashKarp45Adaptive (4, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (4, r, y, err, rr[i] - r);

      BBphic[i]   = y[0];
      BBthetac[i] = y[1];
      qqc[i]      = epsa * rr[i] * BBphi[i] /BBtheta[i];
      PPc[i]      = y[3];

      PPp[i]  = GetPp(rr[i]);
      PPpc[i] = GetPpcrit (rr[i], qqc[i], BBphic[i]);

      if (PPpc[i] > 0.)
	PPpm[i] = PPp[i];
      else
	if (PPp[i] < PPpc[i])
	  PPpm[i] = PPpc[i];
	else
	  PPpm[i] = PPp[i];
     }

  Theta = y[1] /2./y[2];
  Frev  = y[0] /2./y[2];

  printf ("F = %10.4e Theta = %10.4e\n\n", Frev, Theta);

  double Pedge = PPc[Ngrid];
  for (int i = 0; i < Ngrid; i++)
    PPc[i] = PPc[i] - Pedge;
  
  delete[] y; delete[] err;
}

// ############################################
// Function to get equilibrium parallel current
// ############################################
double Pinch::GetSigma (double r)
{
  if (r >= 1.)
    return 0;
  else
    return (2. * epsa /q0) * pow (1. - pow (r, alphas), nus);
}

// #####################################################
// Function to get equilibrium parallel current gradient
// #####################################################
double Pinch::GetSigmap (double r)
{
  if (r >= 1.)
    return 0;
  else
    return - (2. * epsa /q0) * alphas * nus * pow (r, alphas - 1.) * pow (1. - pow (r, alphas), nus - 1.);
}

// ####################################
// Function to get equilibrium pressure
// ####################################
double Pinch::GetP (double r)
{
  if (r >= 1.)
    return 0;
  else
    return (beta0 /2.) * pow (1. - pow (r, alphap), nup);
}

// #############################################
// Function to get equilibrium pressure gradient
// #############################################
double Pinch::GetPp (double r)
{
  if (r >= 1.)
    return 0;
  else
    return - (beta0 /2.) * alphap * nup * pow (r, alphap - 1.) * pow (1. - pow (r, alphap), nup - 1.);
}

// ##############################
// Function to get magnetic shear
// ##############################
double Pinch::Gets (double r, double q)
{
  double sigma = GetSigma (r);

  return 2. - sigma * (epsa*epsa *r*r + q*q) /epsa /q;
}

// ##########################################
// Function to get critical pressure gradient
// ##########################################
double Pinch::GetPpcrit (double r, double q, double Bphi)
{
  double s = Gets (r, q);

  return  - s*s * Bphi*Bphi /r /(1. - q*q);
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Pinch::CashKarp45Rhs (double r, double* y, double* dydr)
{
  if (rhs_chooser == 0)
    {
      double sigma  = GetSigma (r);
      double Pp     = GetPp (r);
      double Bphi   = y[0];
      double Btheta = y[1];
      
      dydr[0] =             - sigma * Btheta - Pp * Bphi   /(Btheta*Btheta + Bphi*Bphi);
      dydr[1] = - Btheta /r + sigma * Bphi   - Pp * Btheta /(Btheta*Btheta + Bphi*Bphi);
      dydr[2] = Bphi * r;
    }
  else
    {
      double sigma  = GetSigma (r);
      double Pp     = GetPp (r);
      double Bphi   = y[0];
      double Btheta = y[1];
      double q      = epsa * r * Bphi /Btheta;
      double Ppc    = GetPpcrit (r, q, Bphi);

      double Ppp;
      if (Ppc > 0.)
	Ppp = Pp;
      else
	if (Pp < Ppc)
	  Ppp = Ppc;
	else
	  Ppp = Pp;

      dydr[0] =             - sigma * Btheta - Ppp * Bphi   /(Btheta*Btheta + Bphi*Bphi);
      dydr[1] = - Btheta /r + sigma * Bphi   - Ppp * Btheta /(Btheta*Btheta + Bphi*Bphi);
      dydr[2] = Bphi * r;
      dydr[3] = Ppp; 
    }
}

