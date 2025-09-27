// Pinch.cpp

#include "Pinch.h"

#define NINPUT 16

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
  bwall  = JSONData["bwall"] .get<double> ();
  dwall  = JSONData["dwall"] .get<double> ();

  mpol = JSONData["mpol"].get<int> ();
  ntor = JSONData["ntor"].get<int> ();
  
  Ngrid = JSONData["Ngrid"].get<int>    ();
  eps   = JSONData["eps"]  .get<double> ();
  delta = JSONData["delta"].get<double> ();

  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  Eta     = JSONData["Eta"]    .get<double> ();
  Maxiter = JSONData["Maxiter"].get<int>    ();
  gmax    = JSONData["gmax"]   .get<double> ();
  Nint    = JSONData["Nint"]   .get<int>    ();

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
  if (bwall < 1.)
    {
      printf ("Pinch:: Error - bwall must be greater than unity\n");
      exit (1);
    }
   if (dwall <= 0.)
    {
      printf ("Pinch:: Error - dwall must be positive\n");
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
  if (delta < 0. || delta > 1.)
    {
      printf ("Pinch:: Error - delta must lie in the interval 0 to 1\n");
      exit (1);
    }
  if (mpol == 0 && ntor == 0)
    {
      printf ("Pinch:: Error - mpol and ntor cannot both be zeron");
      exit (1);
    }
  
  printf ("\n");
  printf ("Class PINCH::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("epsa  = %10.3e q0    = %10.3e beta0 = %10.3e alphas = %10.3e nus = %10.3e alphap = %10.3e nup = %10.3e\n",
	  epsa, q0, beta0, alphas, nus, alphap, nup);
  printf ("bwall = %10.3e dwall = %10.3e mpol  = %3d        ntor   = %3d\n",
	  bwall, dwall, mpol, ntor);
  printf ("Ngrid =  %4d      eps   = %10.3e delta = %10.3e\n",
	  Ngrid, eps, delta);
  
  // ---------------
  // Allocate memory
  // ---------------
  rr       = new double[Ngrid+1];
  rrv      = new double[Ngrid+1];
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

  bbeta1   = new double[Ngrid+1];
  bbeta2   = new double[Ngrid+1];
  ssigp    = new double[Ngrid+1];
  ff       = new double[Ngrid+1];
  gg       = new double[Ngrid+1];

  Bphi_accel   = gsl_interp_accel_alloc ();
  Btheta_accel = gsl_interp_accel_alloc ();
  Pp_accel     = gsl_interp_accel_alloc ();
  q_accel      = gsl_interp_accel_alloc ();
  
  Bphi_spline   = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);
  Btheta_spline = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);
  Pp_spline     = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);
  q_spline      = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);

  Psip   = new double[Ngrid+1];
  Psinw  = new double[Ngrid+1];
  Psipw  = new double[Ngrid+1];
  Psirwm = new double[Ngrid+1];

  // ------------------
  // Set up radial grid
  // ------------------
  double dr = 1. /double (Ngrid);
  for (int i = 0; i <= Ngrid; i++)
    {
      rr[i]  = double (i) * dr;
      rrv[i] = 1. + rr[i];
    }
}

// ###########
// Destructor
// ###########
Pinch::~Pinch ()
{
  delete[] rr;  delete[] ssigma; delete[] PP;     delete[] BBphi;    delete[] BBtheta; delete[] qq;
  delete[] qqc; delete[] PPc;    delete[] BBphic; delete[] BBthetac;
  delete[] PPp; delete[] PPpc;   delete[] PPpm;   delete[] rrv;

  delete[] bbeta1; delete[] bbeta2; delete[] ssigp; delete[] ff; delete[] gg;

  delete[] Psip; delete[] Psinw; delete[] Psipw; delete[] Psirwm; 

  gsl_interp_accel_free (Bphi_accel);  gsl_interp_accel_free (Btheta_accel);  gsl_interp_accel_free (q_accel);
  gsl_interp_accel_free (Pp_accel);
  gsl_spline_free       (Bphi_spline); gsl_spline_free       (Btheta_spline); gsl_spline_free       (q_spline);
  gsl_spline_free       (Pp_spline);
}

// #########################
// Function to solve problem
// #########################
void Pinch::Solve ()
{
  // Set wall data
  epsb = bwall * epsa;
  delw = dwall /bwall;
 
  // Set mode data
  MPOL = double (mpol);
  NTOR = double (ntor);
  ka   = epsa * NTOR;
  kb   = epsb * NTOR;
  qres = MPOL /NTOR;
  
  // Calculate toroidal pinch equilibrium
  CalcEquilibrium ();

  // Find resonant surface
  FindResonant ();

  // Solve Newcomb's equation
  SolveNewcomb ();

  // Calculate vacuum solution
  SolveVacuum ();

  // Calculate resistive wall mode growth-rate
  Growth ();
 
  // Output data to netcdf file
  WriteNetcdf ();
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
  qa    = qq[Ngrid];
  
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

  // ................................
  // Interpolate equilibrium profiles
  // ................................
  gsl_spline_init (Bphi_spline,   rr, BBphic,   Ngrid+1);
  gsl_spline_init (Btheta_spline, rr, BBthetac, Ngrid+1);
  gsl_spline_init (Pp_spline,     rr, PPpm,     Ngrid+1);
  gsl_spline_init (q_spline,      rr, qqc,      Ngrid+1);

  // .............................
  // Calculate additional profiles
  // .............................
  for (int i = 0; i <= Ngrid; i++)
    {
      bbeta1[i] = Getbeta1  (rr[i]);
      bbeta2[i] = Getbeta2  (rr[i]);
      ssigp [i] = GetSigmap (rr[i]);
      ff    [i] = Getf      (rr[i]);
      gg    [i] = tanh(Getg (rr[i]) /10.);
    }
}

// #################################
// Function to find resonant surface
// #################################
void Pinch::FindResonant ()
{
  double x1 = eps;
  double x2 = 1.;
  double F1 = q0 - qres;
  double F2 = qa - qres;
 
  if (F1 * F2 < 0.)
    {
      f_chooser = 0;
      
      Ridder (x1, x2, F1, F2, rres);
    }
  else
    {
      rres = -1.;
    }

  if (rres < 0.)
    printf ("Nonresonant mode\n\n");
  else
    printf ("Rational surface: rres = %10.3e qres = %10.3e residual = %10.3e\n\n", rres, qres, fabs(Getq(rres) - qres));
}

// ####################################
// Function to solve Newcomb's equation
// ####################################
void Pinch::SolveNewcomb ()
{
  for (int i = 0; i <= Ngrid; i++)
    Psip[i] = 0.;

  double  r, h, t_err;
  int     rept; rhs_chooser = 2;
  double* y   = new double[2];
  double* err = new double[2];

  h     = h0;
  count = 0;

  if (rres < 0.)
    {
      r = eps;

      if (mpol == 0)
	{
	  y[0] = 0.;
	  y[1] = 1.;
	}
      else
	{
	  y[0] = pow (eps, fabs (MPOL));
	  y[1] = y[0] /fabs (MPOL);
	}
      
      Psip[0] = 0.;
    }
  else
    {
      r = rres + delta;

      y[0] = 0.;
      y[1] = 1.;
    }

  printf ("Launching solution: r = %10.3e y = (%10.3e, %10.3e)\n", r, y[0], y[1]);
 
  for (int i = 1; i <= Ngrid; i++)
    {
      if (r < rr[i])
	{
	  do
	    {
	      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rr[i]);
	  CashKarp45Fixed (2, r, y, err, rr[i] - r);
	  
	  Psip[i] = y[0];
	}
    }
  lbar = y[1] /Getf (1.) /y[0];
  printf ("Stopping solution: r = %10.3e y = (%10.3e, %10.3e) lbar = %10.3e\n", r, y[0], y[1], lbar);
 
  delete[] y; delete[] err;
  
  double Psia = Psip[Ngrid];
  for (int i = 0; i <= Ngrid; i++)
    Psip[i] /= Psia;
}

// ###########################################
// Function to calculate vacuum eigenfunctions
// ###########################################
void Pinch::SolveVacuum ()
{
  if (ntor == 0)
    {
      Lnw  = - 1./ fabs (MPOL);
      Lpw  = (1. + pow (bwall, - 2. * fabs (MPOL))) / (1. - pow (bwall, - 2. * fabs (MPOL))) /fabs (MPOL);
      alpw = (1. - pow (bwall, - 2. * fabs (MPOL))) /2. /fabs (MPOL);
      
      for (int i = 0; i <= Ngrid; i++)
	{
	  double r = rrv[i];

	  double psinw = pow (r, - fabs (MPOL));
	  
	  double psipw;
	  if (r < bwall)
	    {
	      psipw =   (pow (bwall/r, fabs (MPOL)) - pow (r /bwall, fabs (MPOL)))
		      / (pow (bwall,   fabs (MPOL)) - pow (1./bwall, fabs (MPOL)));
	    }
	  else
	    psipw = 0.;

	  Psinw[i] = psinw;
	  Psipw[i] = psipw;
	}
    }
  else
    {
      double Ima  = GetIm  (1.);
      double Kma  = GetKm  (1.);
      double Impa = GetImp (1.);
      double Kmpa = GetKmp (1.);
      double Impb = GetImp (bwall);
      double Kmpb = GetKmp (bwall);

      Lnw  = - Kma /Kmpa / (fabs (NTOR) * epsa);
      Lpw  = - (Ima * Kmpb - Kma * Impb) /(Impa * Kmpb - Kmpa * Impb) / (fabs (NTOR) * epsa);
      alpw = (Kmpb /Kmpa) * (Impa * Kmpb - Kmpa * Impb) * NTOR*NTOR * epsb*epsb /(MPOL*MPOL + NTOR*NTOR * epsb*epsb);

      for (int i = 0; i <= Ngrid; i++)
	{
	  double r = rrv[i];
	  
	  double Imp = GetImp (r);
	  double Kmp = GetKmp (r);
	  
	  double psinw = r * Kmp /Kmpa;

	  double psipw;
	  if (r < bwall)
	    psipw = r * (Imp * Kmpb - Kmp * Impb) /(Impa * Kmpb - Kmpa * Impb);
	  else
	    psipw = 0.;
	  
	  Psinw[i] = psinw;
	  Psipw[i] = psipw;
	}
    }

  double F = GetF (1.);
  double H = GetH (1.);
  
  Wnw = M_PI*M_PI * F*F * (lbar /H + Lnw);
  Wpw = M_PI*M_PI * F*F * (lbar /H + Lpw);
  c1  =   Wpw / (Wpw - Wnw);
  c2  = - Wnw / (Wpw - Wnw);

  rhs = - Wnw /Wpw /alpw;

  printf ("\nLambda_nw = %10.3e Lambda_pw = %10.3e alpha_w = %10.3e dW_nw = %10.3e dW_pw = %10.3e c1 = %10.3e c2 = %10.3e\n\n",
	  Lnw, Lpw, alpw, Wnw, Wpw, c1, c2);
}

// #####################################################
// Function to calculate resistive wall mode growth-rate
// #####################################################
void Pinch::Growth ()
{
  if (Wpw < 0.)
    {
      gamma = 1.e3;

      printf ("Perfect-wall mode is ideally unstable\n\n"); 
    }
  else
    {
      if (rhs > 0.)
	{
	  f_chooser = 1;
	  
	  gamma = RootFind (0., gmax);
	}
      else
	{
	  f_chooser = 2;
	  
	  gamma = RootFind (0., -gmax);
	}

        printf ("rhs = %10.3e gamma = %10.3e\n\n", rhs, gamma);
    }
}

// ######################################
// Function to output data to netcdf file
// ######################################
void Pinch::WriteNetcdf ()
{
   printf ("Writing data to netcdf file Outputs/Pinch/Pinch.nc:\n");

   double Input[NINPUT];
   
   Input[0]  = epsa;
   Input[1]  = q0;
   Input[2]  = beta0;
   Input[3]  = alphas;
   Input[4]  = nus;
   Input[5]  = alphap;
   Input[6]  = nup;
   Input[7]  = bwall;
   Input[8]  = dwall;
   Input[9]  = MPOL;
   Input[10] = NTOR;
   Input[11] = double (Ngrid);
   Input[12] = eps;
   Input[13] = delta;
   Input[14] = rres;
   Input[15] = qres;
      
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
       NcVar rv_x    = dataFile.addVar ("r_v",             ncDouble, r_d);
       rv_x.putVar (rrv);
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
       NcVar bb1_x   = dataFile.addVar ("beta_1",          ncDouble, r_d);
       bb1_x.putVar (bbeta1);
       NcVar bb2_x   = dataFile.addVar ("beta_2",          ncDouble, r_d);
       bb2_x.putVar (bbeta2);
       NcVar ssp_x   = dataFile.addVar ("sigma_p",         ncDouble, r_d);
       ssp_x.putVar (ssigp);
       NcVar ff_x    = dataFile.addVar ("f",               ncDouble, r_d);
       ff_x.putVar (ff);
       NcVar gg_x    = dataFile.addVar ("g",               ncDouble, r_d);
       gg_x.putVar (gg);
       NcVar psip_x  = dataFile.addVar ("psi_p",           ncDouble, r_d);
       psip_x.putVar (Psip);
       NcVar psinw_x = dataFile.addVar ("psi_nw",          ncDouble, r_d);
       psinw_x.putVar (Psinw);
       NcVar psipw_x = dataFile.addVar ("psi_pw",          ncDouble, r_d);
       psipw_x.putVar (Psipw);
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

// ############################################
// Function to get equilibrium parallel current
// ############################################
double Pinch::GetSigma (double r)
{
  if (r >= 1.)
    return 0.;
  else
    return (2. * epsa /q0) * pow (1. - pow (r, alphas), nus);
}

// #####################################################
// Function to get equilibrium parallel current gradient
// #####################################################
double Pinch::GetSigmap (double r)
{
  if (r >= 1.)
    return 0.;
  else
    return - (2. * epsa /q0) * alphas * nus * pow (r, alphas - 1.) * pow (1. - pow (r, alphas), nus - 1.);
}

// ####################################
// Function to get equilibrium pressure
// ####################################
double Pinch::GetP (double r)
{
  if (r >= 1.)
    return 0.;
  else
    return (beta0 /2.) * pow (1. - pow (r, alphap), nup);
}

// #############################################
// Function to get equilibrium pressure gradient
// #############################################
double Pinch::GetPp (double r)
{
  if (r >= 1.)
    return 0.;
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

  return - s*s * Bphi*Bphi /r /(1. - q*q) /8.;
}

// #################
// Function to get q
// #################
double Pinch::Getq (double r)
{
  if (r < 1.)
    {
      return gsl_spline_eval (q_spline, r, q_accel);
    }
  else
    {
      return qq[Ngrid] * r*r;
    }
}

// ######################
// Function to get Btheta
// ######################
double Pinch::GetBtheta (double r)
{
  if (r < 1.)
    {
      return gsl_spline_eval (Btheta_spline, r, Btheta_accel);
    }
  else
    {
      return BBtheta[Ngrid] /r;
    }
}

// ####################
// Function to get Bphi
// ####################
double Pinch::GetBphi (double r)
{
  if (r < 1.)
    {
      return gsl_spline_eval (Bphi_spline, r, Bphi_accel);
    }
  else
    {
      return BBphi[Ngrid];
    }
}

// #####################
// Function to get beta1
// #####################
double Pinch::Getbeta1 (double r)
{
  if (r < 1.)
    {
      double Pp = gsl_spline_eval (Pp_spline, r, Pp_accel);

      return r * Pp;
    }
  else
    return 0.;
 }

// #####################
// Function to get beta2
// #####################
double Pinch::Getbeta2 (double r)
{
  if (r < 1.)
    {
      double Pp2 = gsl_spline_eval_deriv (Pp_spline, r, Pp_accel);

      return r*r * Pp2;
    }
  else
    return 0.;
}   

// #################
// Function to get F
// #################
double Pinch::GetF (double r)
{
  double Btheta = GetBtheta (r);
  double Bphi   = GetBphi   (r);
  
  return MPOL * Btheta - NTOR * epsa * r * Bphi;
}

// #################
// Function to get G
// #################
double Pinch::GetG (double r)
{
  double Btheta = GetBtheta (r);
  double Bphi   = GetBphi   (r);

  return MPOL * Bphi + NTOR * epsa * r * Btheta;
}

// #################
// Function to get H
// #################
double Pinch::GetH (double r)
{
  return MPOL*MPOL + NTOR*NTOR * epsa*epsa * r*r;
}

// #################
// Function to get f
// #################
double Pinch::Getf (double r)
{
  double H = GetH (r);

  return 1. /H;
}

// #################
// Function to get g
// #################
double Pinch::Getg (double r)
{
  double Bt  = GetBtheta (r);
  double Bp  = GetBphi   (r);
  double sg  = GetSigma  (r);
  double sgp = GetSigmap (r);
  double b1  = Getbeta1  (r);
  double b2  = Getbeta2  (r);
  double q   = Getq      (r);
  double F   = GetF      (r);
  double G   = GetG      (r);
  double H   = GetH      (r);

  double B2   = Bt*Bt + Bp*Bp;
  double B4   = B2*B2;
  double m    = MPOL;
  double neps = NTOR * epsa * r;

  return
    1.
    + 2.*m*neps*r*sg /H/H
    - (r*r * sg*sg  + (b2 - b1) /B2 + b1 * (b1 + 2.*Bt*Bt) /B4) /H
    + (r*r * sgp - 2.*r*sg*b1 /B2 + 2.*m*neps*b1 /H/B2 + 2.*neps*Bt*b1*(1. - q*q) /F/B2) * G /H/F;
}

// ##################
// Function to get Im
// ##################
double Pinch::GetIm (double r)
{
  return gsl_sf_bessel_In (abs (mpol), r * ka);
}

// ##################
// Function to get Km
// ##################
double Pinch::GetKm (double r)
{
  return gsl_sf_bessel_Kn (abs (mpol), r * ka);
}

// ###################
// Function to get Imp
// ###################
double Pinch::GetImp (double r)
{
  return gsl_sf_bessel_In (abs (mpol) + 1, r * ka) + (abs (mpol) /(r * ka)) * gsl_sf_bessel_In (abs (mpol), r * ka);
}

// ###################
// Function to get Kmp
// ###################
double Pinch::GetKmp (double r)
{
  return - gsl_sf_bessel_Kn (abs (mpol) + 1, r * ka) + (abs (mpol) /(r * ka)) * gsl_sf_bessel_Kn (abs (mpol), r * ka);
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
  else if (rhs_chooser == 1)
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
  else
    {
      double f = Getf (r);
      double g = Getg (r);

      dydr[0] = y[1] /r/f;
      dydr[1] = g * y[0] /r;
    }
}

// ################################
// Target function for zero finding
// ################################
double Pinch::RootFindF (double x)
{
  if (f_chooser == 0)
    {
      double q = Getq (x);

      return q - qres;
    }
  else if (f_chooser == 1)
    {
      return sqrt (x /delw)   * tanh (sqrt (delw * x))  - rhs;
    }
  else if (f_chooser == 2)
    {
      return sqrt (- x /delw) * tan (sqrt (- delw * x)) + rhs;
    }
}
