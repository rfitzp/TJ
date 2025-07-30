// ECE.cpp

#include "ECE.h"

#define NINPUT 6

// ###########
// Constructor
// ###########
ECE::ECE ()
{
  // ------------------------------------------
  // Ensure that directory ../Outputs/ECE exits
  // ------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/ECE"))
    {
      exit (1);
    }

  // ---------------------------------------
  // Read control parameters from JSON files
  // ---------------------------------------
  string JSONFilename = "../Inputs/ECE.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  Te = JSONData["Te"].get<double>();
  ne = JSONData["ne"].get<double>();
  B0 = JSONData["B0"].get<double>();
  R0 = JSONData["R0"].get<double>();
  Rw = JSONData["Rw"].get<double>();

  zmax  = JSONData["zmax"] .get<double>();
  znum  = JSONData["znum"] .get<int>();
  wcmin = JSONData["wcmin"].get<double>();
  wcmax = JSONData["wcmax"].get<double>();
  wcnum = JSONData["wcnum"].get<int>();

  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (Te < 0.)
    {
      printf ("ECE: Error - Te must be positive\n");
      exit (1);
    }
  if (ne < 0.)
    {
      printf ("ECE: Error - ne must be positive\n");
      exit (1);
    }
  if (B0 < 0.)
    {
      printf ("ECE: Error - B0 must be positive\n");
      exit (1);
    }
  if (R0 < 0.)
    {
      printf ("ECE: Error - R0 must be positive\n");
      exit (1);
    }
  if (Rw < 0.)
    {
      printf ("ECE: Error - Rw must be positive\n");
      exit (1);
    }

  // ---------------------
  // Print welcome message
  // ---------------------
  printf ("\nClass ECE::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
}

// ##########
// Destructor
// ##########
ECE::~ECE ()
{
}

// #########################
// Function to solve problem
// #########################
void ECE::Solve ()
{
  // ----------------------------
  // Calculate derived parameters
  // ----------------------------
  SetDerivedParameters (Te, ne, B0, R0, Rw);
  
  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("Calculation Parameters:\n");
  printf ("Te  = %11.4e ne   = %11.4e B0   = %11.4e R0   = %11.4e Rw  = %11.4e\n",
	  Te, ne, B0, R0, Rw);
  printf ("wc0 = %11.4e tau0 = %11.4e vtc2 = %11.4e w12  = %11.4e w22 = %11.4e tw = %11.4e\n",
	  wc0, tau0, vtc2, w1*w1, 4.*w1*w1, tw);
  printf ("acc = %11.4e h0   = %11.4e hmin = %11.4e hmax = %11.4e\n",
	  acc, h0, hmin, hmax);

  // -----------------------------------------
  // Calculate F_(7/2)(z) on array of z values
  // -----------------------------------------
  zz   = new double [znum];
  F72r = new double [znum];
  F72i = new double [znum];
  
  for (int i = 0; i < znum; i++)
    {
      double          z   = - zmax + 2.*zmax * double (i) /double (znum-1);
      complex<double> F72 = Get_F72 (z);
      
      zz  [i] = z;
      F72r[i] = real (F72);
      F72i[i] = imag (F72);
    }

  // ---------------------------------
  // Calculate absorption coefficients
  // ---------------------------------
  wwc    = new double [wcnum];
  alphaO = new double [wcnum];
  alphaX = new double [wcnum];

  for (int i = 0; i < wcnum; i++)
    {
      double wc  = wcmin * w1 + (wcmax - wcmin) * w1 * double (i) /double (wcnum-1);
      double alO = Get_alpha1O (wc);
      double alX = Get_alpha2X (wc);

      wwc   [i] = wc;
      alphaO[i] = alO;
      alphaX[i] = alX;
    }

  // ------------------------
  // Calculate optical depths
  // ------------------------
  tauO = new double [wcnum];
  tauX = new double [wcnum];

  double wc, h, t_err;
  int    rept; count = 0; rhs_chooser = 0;

  double* y     = new double[2];
  double* dydwc = new double[2];
  double* err   = new double[2];
 
  for (int i = 0; i < wcnum; i++)
    {
      h     = h0;
      count = 0;
      wc    = w1;
      y[0]  = 0.;
      y[1]  = 0.;
      
      do
	{
	  CashKarp45Adaptive (2, wc, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (wc < wwc[i]);
      CashKarp45Fixed (2, wc, y, err, wwc[i] - wc);

      tauO[i] = y[0];
      tauX[i] = y[1];
    }

  tauOst = tauO[wcnum-1];
  tauXst = tauX[wcnum-1];

  delete[] y; delete[] dydwc; delete[] err;

  // ----------------------------------------
  // Calculate spectral convolution functions
  // ----------------------------------------
  HO = new double[wcnum];
  HX = new double[wcnum];
  RR = new double[wcnum];
  FO = new double[wcnum];
  FX = new double[wcnum];
  
  for (int i = 0; i < wcnum; i++)
    {
      HO[i] = tau0 * alphaO[i] * exp(- tauO[i]) /wwc[i] /(1. - exp (- tauOst));
      HX[i] = tau0 * alphaX[i] * exp(- tauX[i]) /wwc[i] /(1. - exp (- tauXst));

      RR[i] = R0 * wc0 /wp/ wwc[i];
      FO[i] = (wc0 /wp) * HO[i] * R0 /RR[i]/RR[i];
      FX[i] = (wc0 /wp) * HX[i] * R0 /RR[i]/RR[i];
    }

  // --------------------------------------------------------------------------------------------------------
  // Determine standard deviations and radial shifts of spatial convolution functions via cubic interpolation
  // --------------------------------------------------------------------------------------------------------
  int    iO    = 0,     iX    = 0;
  double FmaxO = -1.e6, FmaxX = -1.e6;

  for (int i = 0; i < wcnum; i++)
    {
      if (FO[i] > FmaxO)
	{
	  FmaxO = FO[i];
	  iO    = i;
	}
      if (FX[i] > FmaxX)
	{
	  FmaxX = FX[i];
	  iX    = i;
	}
    }

  if (iO == 0)
    iO = 1;
  else if (iO == wcnum - 1)
    iO = wcnum - 2;
  if (iX == 0)
    iX = 1;
  else if (iX == wcnum - 1)
    iX = wcnum - 2;

  double xm1, x0, x1, ym1, y0, y1, fm1, f0, f1, num, denom, x, max;

  xm1 = RR[iO-1];
  x0  = RR[iO  ];
  x1  = RR[iO+1];
  ym1 = FO[iO-1];
  y0  = FO[iO  ];
  y1  = FO[iO+1];

  fm1 = (xm1 - x0)  * (xm1 - x1);
  f0  = (x0  - xm1) * (x0  - x1);
  f1  = (x1  - xm1) * (x1  - x0);

  num   = ym1 * (x0 + x1) /fm1 + y0 * (xm1 + x1) /f0 + y1 * (xm1 + x0) /f1;
  denom = 2. * (ym1 /fm1 + y0 /f0 + y1 /f1);

  x      = num /denom;
  DeltaO = Rw - x;
  max    = ym1 * (x - x0) * (x - x1) /fm1 + y0 * (x - xm1) * (x - x1) /f0 + y1 * (x - xm1) * (x - x0) /f1;
  sigmaO = 1. /mc_sp2 /max;

  xm1 = RR[iX-1];
  x0  = RR[iX  ];
  x1  = RR[iX+1];
  ym1 = FX[iX-1];
  y0  = FX[iX  ];
  y1  = FX[iX+1];

  fm1 = (xm1 - x0)  * (xm1 - x1);
  f0  = (x0  - xm1) * (x0  - x1);
  f1  = (x1  - xm1) * (x1  - x0);

  num   = ym1 * (x0 + x1) /fm1 + y0 * (xm1 + x1) /f0 + y1 * (xm1 + x0) /f1;
  denom = 2. * (ym1 /fm1 + y0 /f0 + y1 /f1);

  x      = num /denom;
  DeltaX = Rw - x;
  max    = ym1 * (x - x0) * (x - x1) /fm1 + y0 * (x - xm1) * (x - x1) /f0 + y1 * (x - xm1) * (x - x0) /f1;
  sigmaX = 1. /mc_sp2 /max;

  printf ("\nResults:\n");
  printf ("Delta_1^O = %11.4e m sigma_1^O = %11.4e m tau_1^O = %11.4e\n", DeltaO, sigmaO, tauOst);
  printf ("Delta_2^X = %11.4e m sigma_2^X = %11.4e m tau_2^X = %11.4e\n", DeltaX, sigmaX, tauXst);

  // ------------------------------------------------------------------
  // Calculate truncated Gaussian fits to spatial convolution functions
  // ------------------------------------------------------------------
  FOfit = new double[wcnum];
  FXfit = new double[wcnum];

  for (int i = 0; i < wcnum; i++)
    {
      double R   = RR[i];
      double xiO = (R - Rw + DeltaO) /sigmaO;
      double xiX = (R - Rw + DeltaX) /sigmaX;

      double PO = gsl_cdf_gaussian_P (DeltaO /sigmaO, 1.);
      double PX = gsl_cdf_gaussian_P (DeltaX /sigmaX, 1.);

      if (R > Rw)
	{
	  FOfit[i] = 0.;
	  FXfit[i] = 0.;
	}
      else
	{
	  FOfit[i] = exp ( - xiO*xiO /2.) /mc_sp2 /sigmaO /PO;
	  FXfit[i] = exp ( - xiX*xiX /2.) /mc_sp2 /sigmaX /PX;
	}
    }

  // ------------------------------------------
  // Calculate Gaussian fit parameters directly
  // ------------------------------------------
  double sigma, Delta, taust;
  printf ("\nDirect Calculation:\n");

  GetOmodeFit (Te, ne, B0, R0, Rw, sigma, Delta, taust);
  printf ("Delta_1^O = %11.4e m sigma_1^O = %11.4e m tau_1^0 = %11.4e\n",   Delta, sigma, taust);

  GetXmodeFit (Te, ne, B0, R0, Rw, sigma, Delta, taust);
  printf ("Delta_2^X = %11.4e m sigma_2^X = %11.4e m tau_2^X = %11.4e\n\n", Delta, sigma, taust);

  // -----------------
  // Write netcdf file
  // -----------------
  Write_netcdf ();

  // --------
  // Clean up
  // --------
  delete[] zz;     delete[] F72r; delete[] F72i;
  delete[] wwc;    delete[] RR;
  delete[] alphaO; delete[] tauO; delete[] HO; delete[] FO; delete[] FOfit;
  delete[] alphaX; delete[] tauX; delete[] HX; delete[] FX; delete[] FXfit;
}

// ##################################
// Function to set derived parameters
// ##################################
void ECE::SetDerivedParameters (double _Te, double _ne, double _B0, double _R0, double _Rw)
{
  wc0  = pc_e * _B0 /pc_m_e;
  tau0 = wc0 * _R0 /pc_c;
  vt   = sqrt (_Te * pc_e /pc_m_e);
  vtc2 = vt*vt /pc_c/pc_c;
  wp   = sqrt (_ne * pc_e*pc_e /pc_epsilon_0 /pc_m_e);
  w1   = (wc0 /wp) * (_R0 /_Rw);
  tw   = _Te * pc_e /pc_m_e /pc_c/pc_c;
}

// ##############################################################################################################
// Function to calculate standard deviation and radial shift of 1st harmonic O-mode spatial convolution functions
// ##############################################################################################################
void ECE::GetOmodeFit (double _Te, double _ne, double _B0, double _R0, double _Rw, double& sigma, double& Delta, double& taust)
{
  // ----------------------
  // Set derived parameters
  // ----------------------
  SetDerivedParameters (_Te, _ne, _B0, _R0, _Rw);

  // ------------------------------------------
  // Find peak of spatial distribution function
  // ------------------------------------------
  double wc, h, t_err;
  int    rept; count = 0; rhs_chooser = 1;

  double* y     = new double[1];
  double* dydwc = new double[1];

  double R = _Rw, G, F = 0., xm1, x0 = _Rw, x1, ym1, y0 = 0., y1, fm1, f0, f1, num, denom, x, max;

  h     = h0;
  count = 0;
  wc    = w1;
  y[0]  = 0.;
  
  do
    {
      xm1 = x0;
      ym1 = y0;
      x0  = R;
      y0  = F;
 
      CashKarp45Adaptive (2, wc, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);

      R = _R0 * wc0 /wp /wc;
      G = tau0 * Get_alpha1O (wc) * exp (- y[0]) /wc;
      F = (wc0 /wp) * G * _R0 /R/R;

      x1 = R;
      y1 = F;
    }
  while (y1 > y0);

  // ---------------
  // Calculate Delta
  // ---------------
  fm1 = (xm1 - x0)  * (xm1 - x1);
  f0  = (x0  - xm1) * (x0  - x1);
  f1  = (x1  - xm1) * (x1  - x0);

  num   = ym1 * (x0 + x1) /fm1 + y0 * (xm1 + x1) /f0 + y1 * (xm1 + x0) /f1;
  denom = 2. * (ym1 /fm1 + y0 /f0 + y1 /f1);

  x     = num /denom;
  Delta = _Rw - x;
  max   = ym1 * (x - x0) * (x - x1) /fm1 + y0 * (x - xm1) * (x - x1) /f0 + y1 * (x - xm1) * (x - x0) /f1;

  // ---------------
  // Calculate taust
  // ---------------
  do
    {
      CashKarp45Adaptive (2, wc, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (wc < wcmax * w1);

  taust = y[0];
  
  // ---------------
  // Calculate sigma
  // ---------------
  sigma = (1. - exp (- taust)) /mc_sp2 /max;
  
  delete[] y; delete[] dydwc;
}

// ##############################################################################################################
// Function to calculate standard deviation and radial shift of 2nd harmonic X-mode spatial convolution functions
// ##############################################################################################################
void ECE::GetXmodeFit (double _Te, double _ne, double _B0, double _R0, double _Rw, double& sigma, double& Delta, double& taust)
{
  // ----------------------
  // Set derived parameters
  // ----------------------
  SetDerivedParameters (_Te, _ne, _B0, _R0, _Rw);

  // ------------------------------------------
  // Find peak of spatial distribution function
  // ------------------------------------------
  double wc, h, t_err;
  int    rept; count = 0; rhs_chooser = 2;

  double* y     = new double[1];
  double* dydwc = new double[1];

  double R = _Rw, G, F = 0., xm1, x0 = _Rw, x1, ym1, y0 = 0., y1, fm1, f0, f1, num, denom, x, max;

  h     = h0;
  count = 0;
  wc    = w1;
  y[0]  = 0.;
  
  do
    {
      xm1 = x0;
      ym1 = y0;
      x0  = R;
      y0  = F;
 
      CashKarp45Adaptive (2, wc, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);

      R = _R0 * wc0 /wp /wc;
      G = tau0 * Get_alpha2X (wc) * exp (- y[0]) /wc;
      F = (wc0 /wp) * G * _R0 /R/R;

      x1 = R;
      y1 = F;
    }
  while (y1 > y0);

  // ---------------
  // Calculate Delta
  // ---------------
  fm1 = (xm1 - x0)  * (xm1 - x1);
  f0  = (x0  - xm1) * (x0  - x1);
  f1  = (x1  - xm1) * (x1  - x0);

  num   = ym1 * (x0 + x1) /fm1 + y0 * (xm1 + x1) /f0 + y1 * (xm1 + x0) /f1;
  denom = 2. * (ym1 /fm1 + y0 /f0 + y1 /f1);

  x     = num /denom;
  Delta = _Rw - x;
  max   = ym1 * (x - x0) * (x - x1) /fm1 + y0 * (x - xm1) * (x - x1) /f0 + y1 * (x - xm1) * (x - x0) /f1;

  // ---------------
  // Calculate taust
  // ---------------
  do
    {
      CashKarp45Adaptive (2, wc, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (wc < wcmax * w1);

  taust = y[0];

  // ---------------
  // Calculate sigma
  // ---------------
  sigma = (1. - exp (- taust)) /mc_sp2 /max;
  
  delete[] y; delete[] dydwc;
}

// ########################
// Function to evaluate F72
// ########################
complex<double> ECE::Get_F72 (double z)
{
  complex<double> sqz;
  
  if (z > 0)
    sqz = sqrt (z);
  else
    sqz = sqrt (-z) * mc_i;
        
  complex<double> f72 = (8./15.) * (0.75 - 0.5*z + z*z + mc_i * z*z*sqz * ZPlasma (mc_i * sqz));

  return f72;
}

// ######################
// Get resonance function
// ######################
double ECE::Get_z (double wc)
{
  return (1. - wc/w1) /vtc2;
}

// ###########################################
// Get refractive index of 1st harmonic O-mode
// ###########################################
double ECE::Get_NperpO (double wc)
{
  double          z    = Get_z   (wc);
  complex<double> F72  = Get_F72 (z);
  double          F72r = real    (F72);

  double wpc2  = 1. /wc/wc;
  double denom = 1. + 0.5 * wpc2 * F72r;
  double N2    = (1. - wpc2) /denom;

  if (N2 > 0.)
    return sqrt (N2);
  else
    return 0.;
}

// ###########################################
// Get refractive index of 2nd harmonic X-mode
// ###########################################
double ECE::Get_NperpX (double wc)
{
  double          z    = Get_z   (wc);
  complex<double> F72  = Get_F72 (z);
  double          F72r = real    (F72);

  double wpc2 = 1. /wc/wc;
  double NC2  = 1. - (wpc2 /3.) * (1. - 0.25 * wpc2) /(1. - wpc2 /3.);
  double a    = - 0.5 * wpc2 * F72r / (1. - wpc2 /3.);
  double b    = - 2. * (1. - wpc2 /6.) * a;
  double N2   = NC2 * (1. - (b + a * NC2));
    
  if (N2 > 0.)
    return sqrt (N2);
  else
    return 0.;
}

// #################################################
// Get Absorption coefficient of 1st harmonic O-mode
// #################################################
double ECE::Get_alpha1O (double wc)
{
  double          z    = Get_z   (wc);
  complex<double> F72  = Get_F72 (z);
  double          F72r = real    (F72);
  double          F72i = imag    (F72);

  double wpc2   = 1. /wc/wc;
  double denom  = 1. + 0.5 * wpc2 * F72r;
  double NperpO = Get_NperpO (wc);
 
  return 0.5 * NperpO * wpc2 * (-F72i) /denom;
}

// #################################################
// Get absorption coefficient of 2nd harmonic X-mode
// #################################################
double ECE::Get_alpha2X (double wc)
{
  double          z    = Get_z   (wc);
  complex<double> F72  = Get_F72 (z);
  double          F72r = real    (F72);
  double          F72i = imag    (F72);
 
  double  wpc2   = 1. /wc/wc;
  double  NperpX = Get_NperpX (wc);
  double  N2     = NperpX * NperpX;
  double  a2     = 0.5 * wpc2 * (1. + 3. * N2 * F72r) /(3. - wpc2 * (1. + 1.5 * N2 * F72r));
  double  A2     = (1. + a2) * (1. + a2);
  double  denom  = 1. + 0.5 * wpc2 * A2 * F72r;

  return NperpX * wpc2 * A2 * (-F72i) /denom;
}

// ###################################
// Right hand sides of layer equations
// ###################################
void ECE::CashKarp45Rhs (double wc, double* y, double* dydwc)
{
  if (rhs_chooser == 0)
    {
      if (w1 < wc)
	{
	  dydwc[0]= tau0 * Get_alpha1O (wc) /wc;
	  dydwc[1]= tau0 * Get_alpha2X (wc) /wc;
	}
      else
	{
	  dydwc[0] = 0.;
	  dydwc[1] = 0.;
	}
    }
  else if (rhs_chooser == 1)
    {
      if (w1 < wc)
	{
	  dydwc[0]= tau0 * Get_alpha1O (wc) /wc;
	}
      else
	{
	  dydwc[0] = 0.;
	}
    }
  else
    {
      if (w1 < wc)
	{
	  dydwc[0]= tau0 * Get_alpha2X (wc) /wc;
	}
      else
	{
	  dydwc[0] = 0.;
	}
    }
}

// #############################
// Function to write netcdf file
// #############################
void ECE::Write_netcdf ()
{
  printf ("Writing data to netcdf file Outputs/ECE/ECE.nc:\n");

   double Input[NINPUT];
   
   Input[0] = Te;
   Input[1] = ne;
   Input[2] = B0;
   Input[3] = R0;
   Input[4] = Rw;
   Input[5] = w1;
   
   try
     {
       NcFile dataFile ("../Outputs/ECE/ECE.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim i_d = dataFile.addDim ("Ni", NINPUT);
       NcDim z_d = dataFile.addDim ("Nz", znum);
       NcDim w_d = dataFile.addDim ("Nc", wcnum);

       NcVar i_x    = dataFile.addVar ("InputParameters", ncDouble, i_d);
       i_x.putVar (Input);
       NcVar z_x    = dataFile.addVar ("z",               ncDouble, z_d);
       z_x.putVar (zz);
       NcVar F72r_x = dataFile.addVar ("F72_r",           ncDouble, z_d);
       F72r_x.putVar (F72r);
       NcVar F72i_x = dataFile.addVar ("F72_i",           ncDouble, z_d);
       F72i_x.putVar (F72i);
       NcVar wc_x   = dataFile.addVar ("wc",              ncDouble, w_d);
       wc_x.putVar (wwc);
       NcVar R_x    = dataFile.addVar ("R",               ncDouble, w_d);
       R_x.putVar (RR);
       NcVar alO_x  = dataFile.addVar ("alpha_1^O",       ncDouble, w_d);
       alO_x.putVar (alphaO);
       NcVar alX_x  = dataFile.addVar ("alpha_2^X",       ncDouble, w_d);
       alX_x.putVar (alphaX);
       NcVar tauO_x = dataFile.addVar ("tau_1^O",         ncDouble, w_d);
       tauO_x.putVar (tauO);
       NcVar tauX_x = dataFile.addVar ("tau_2^X",         ncDouble, w_d);
       tauX_x.putVar (tauX);
       NcVar HO_x   = dataFile.addVar ("H_1^O",           ncDouble, w_d);
       HO_x.putVar (HO);
       NcVar HX_x   = dataFile.addVar ("H_2^X",           ncDouble, w_d);
       HX_x.putVar (HX);
       NcVar FO_x   = dataFile.addVar ("F_1^O",           ncDouble, w_d);
       FO_x.putVar (FO);
       NcVar FX_x   = dataFile.addVar ("F_2^X",           ncDouble, w_d);
       FX_x.putVar (FX);
       NcVar FOf_x  = dataFile.addVar ("F_1^O_fit",       ncDouble, w_d);
       FOf_x.putVar (FOfit);
       NcVar FXf_x  = dataFile.addVar ("F_2^X_fit",       ncDouble, w_d);
       FXf_x.putVar (FXfit);
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/ECE/ECE.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }
}
