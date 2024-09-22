// TJ.cpp

#include "TJ.h"

// ###########
// Constructor
// ###########
TJ::TJ ()
{
  // ---------------------------
  // Set root finding parameters
  // ---------------------------
  Eta     = 1.e-16;
  Maxiter = 30;
  
  // -----------------------------------
  // Set adaptive integration parameters
  // -----------------------------------
  maxrept = 50;
  flag    = 2;
  
  // --------------------------------
  // Set Cash-Karp RK4/RK5 parameters
  // --------------------------------
  aa1  = 0.;
  aa2  = 1./5.;
  aa3  = 3./10.;
  aa4  = 3./5.;
  aa5  = 1.;
  aa6  = 7./8.;

  cc1  =  37./378.;
  cc3  = 250./621.;
  cc4  = 125./594.;
  cc6  = 512./1771.;

  ca1  = cc1 -  2825./27648.;
  ca3  = cc3 - 18575./48384.;
  ca4  = cc4 - 13525./55296.;
  ca5  =     -   277./14336.;
  ca6  = cc6 -     1./4.;

  bb21 = 1./5.;

  bb31 = 3./40.;
  bb32 = 9./40.;

  bb41 =   3./10.;
  bb42 = - 9./10.;
  bb43 =   6./5.;

  bb51 = - 11./54.;
  bb52 =    5./2.;
  bb53 = - 70./27.;
  bb54 =   35./27.;

  bb61 =  1631./55296.;
  bb62 =   175./512.;
  bb63 =   575./13824.;
  bb64 = 44275./110592.;
  bb65 =   253./4096.;

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "Inputs/TJ.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR    = JSONData["TJ_control"]["NTOR"]   .get<int> ();
  MMIN    = JSONData["TJ_control"]["MMIN"]   .get<int> ();
  MMAX    = JSONData["TJ_control"]["MMAX"]   .get<int> ();
  EPS     = JSONData["TJ_control"]["EPS"]    .get<double> ();
  DEL     = JSONData["TJ_control"]["DEL"]    .get<double> ();
  NFIX    = JSONData["TJ_control"]["NFIX"]   .get<int> ();
  NDIAG   = JSONData["TJ_control"]["NDIAG"]  .get<int> ();
  NULC    = JSONData["TJ_control"]["NULC"]   .get<double> ();
  ITERMAX = JSONData["TJ_control"]["ITERMAX"].get<int> ();
  FREE    = JSONData["TJ_control"]["FREE"]   .get<int> ();
  acc     = JSONData["TJ_control"]["acc"]    .get<double> ();
  h0      = JSONData["TJ_control"]["h0"]     .get<double> ();
  hmin    = JSONData["TJ_control"]["hmin"]   .get<double> ();
  hmax    = JSONData["TJ_control"]["hmax"]   .get<double> ();
  EPSF    = JSONData["TJ_control"]["EPSF"]   .get<double> ();
  B0      = JSONData["TJ_control"]["B0"]     .get<double> ();
  R0      = JSONData["TJ_control"]["R0"]     .get<double> ();
  n0      = JSONData["TJ_control"]["n0"]     .get<double> ();
  alpha   = JSONData["TJ_control"]["alpha"]  .get<double> ();
  Zeff    = JSONData["TJ_control"]["Zeff"]   .get<double> ();
  Mion    = JSONData["TJ_control"]["Mion"]   .get<double> ();
  Chip    = JSONData["TJ_control"]["Chip"]   .get<double> ();
  Teped   = JSONData["TJ_control"]["Teped"]  .get<double> ();

  // ............
  // Sanity check
  // ............
  if (NTOR < 1)
    {
      printf ("TJ: Error - NTOR must be positive\n");
      exit (1);
    }
  if (MMAX < MMIN)
    {
      printf ("TJ: Error - MMIN must be less that MMAX\n");
      exit (1);
    }
  if (MMAX < MMIN)
    {
      printf ("TJ: Error - MMIN must be less that MMAX\n");
      exit (1);
    }
  if (EPS <= 0.)
    {
      printf ("TJ: Error - EPS must be positive\n");
      exit (1);
    }
  if (DEL <= 0.)
    {
      printf ("TJ: Error - DEL must be positive\n");
      exit (1);
    }
  if (NFIX < 0)
    {
      printf ("TJ: Error - NFIX cannot be negative\n");
      exit (1);
    }
  if (NDIAG < 0)
    {
      printf ("TJ: Error - NDIAG cannot be less that two\n");
      exit (1);
    }
  if (NULC <= 0.)
    {
      printf ("TJ: Error - NULC must be positive\n");
      exit (1);
    }
  if (ITERMAX < 0)
    {
      printf ("TJ: Error - ITERMAX cannot be negative\n");
      exit (1);
    }
    if (acc <= 0.)
    {
      printf ("TJ:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("TJ:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("TJ:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("TJ:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("TJ:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (EPSF <= 0.)
    {
      printf ("TJ: Error - EPSF must be positive\n");
      exit (1);
    }
  if (R0 <= 0.)
    {
      printf ("TJ: Error - R0 must be positive\n");
      exit (1);
    }
  if (n0 <= 0.)
    {
      printf ("TJ: Error - n0 must be positive\n");
      exit (1);
    }
  if (Zeff < 1.)
    {
      printf ("TJ: Error - Zeff cannot be less than unity\n");
      exit (1);
    }
  if (Mion < 1.)
    {
      printf ("TJ: Error - Mion cannot be less than unity\n");
      exit (1);
    }
  if (Chip <= 0.)
    {
      printf ("TJ: Error - Chip must be positive\n");
      exit (1);
    }
  if (Teped <= 0.)
    {
      printf ("TJ: Error - Teped must be positive\n");
      exit (1);
    }
  
  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\nClass TJ::\n");
  printf ("Calculation parameters:\n");
  printf ("ntor = %3d        mmin  = %3d        mmax = %3d        eps     = %10.3e del  = %10.3e\n",
	  NTOR, MMIN, MMAX, EPS, DEL);
  printf ("nfix = %3d        ndiag = %3d       nulc = %10.3e itermax = %3d        free =  %1d\n",
	  NFIX, NDIAG, NULC, ITERMAX, FREE);
  printf ("acc  = %10.3e h0    = %10.3e hmin = %10.3e hmax    = %10.3e epsf = %10.3e\n",
	  acc, h0, hmin, hmax, EPSF);
  printf ("B0   = %10.3e R0    = %10.3e n0   = %10.3e alpha   = %10.3e Zeff = %10.3e Mion = %10.3e Chip = %10.3e Teped = %10.3e\n",
	  B0, R0, n0, alpha, Zeff, Mion, Chip, Teped);
}

// ##########
// Destructor
// ##########
TJ::~TJ ()
{
}

// #########################
// Function to solve problem
// #########################
void TJ::Solve ()
{ 
  // Set toroidal and poloidal mode numbers
  SetModeNumbers ();

  // Read equilibrium data
  ReadEquilibrium ();

  // Read RMP coil data
  ReadCoils ();

  // Calculate metric data at plasma boundary
  CalculateMetric ();

  // Calculate vacuum matrices
  GetVacuum ();

  // Find rational surfaces
  FindRational ();

  // Calculate resonant layer data
  GetLayerData ();

  // Solve outer region odes
  ODESolve ();

  // Determine tearing mode dispersion relation and tearing eigenfunctions
  FindDispersion ();

  // Calculate resonant magnetic perturbation data
  CalculateResonantMagneticPerturbation ();

  // Calculate unreconnected eigenfunction and RMP response visualization data
  VisualizeEigenfunctions ();

  // Calculate ideal stability
  CalculateIdealStability ();
  
  // Write program data to Netcdf file
  WriteNetcdf ();

  // Clean up
  CleanUp ();
}

// ##################################################
// Function to set toroidal and poloidal mode numbers
// ##################################################
void TJ::SetModeNumbers ()
{
  ntor = double (NTOR);
  J    = MMAX - MMIN + 1;
  MPOL = new int[J];
  mpol = new double[J];

  for (int j = 0; j < J; j++)
    {
      MPOL[j] = MMIN + j;
      mpol[j] = double (MMIN + j);
    }
}

// #############################
// Function to deallocate memory
// #############################
void TJ::CleanUp ()
{
  printf ("Cleaning up:\n");
  
  delete[] rr; delete[] pp; delete[] ppp; delete[] q; 
  delete[] s;  delete[] s2; delete[] S1;  delete[] P1;
  delete[] P2; delete[] P3; delete[] g2;  delete[] p2;

  gsl_spline_free (g2spline);
  gsl_spline_free (p2spline);
  gsl_spline_free (ppspline);
  gsl_spline_free (pppspline);
  gsl_spline_free (qspline);
  gsl_spline_free (sspline);
  gsl_spline_free (s2spline);
  gsl_spline_free (S1spline);
  gsl_spline_free (P1spline);
  gsl_spline_free (P2spline);
  gsl_spline_free (P3spline);

  gsl_interp_accel_free (g2acc);
  gsl_interp_accel_free (p2acc);
  gsl_interp_accel_free (ppacc);
  gsl_interp_accel_free (pppacc);
  gsl_interp_accel_free (qacc);
  gsl_interp_accel_free (sacc);
  gsl_interp_accel_free (s2acc);
  gsl_interp_accel_free (S1acc);
  gsl_interp_accel_free (P1acc);
  gsl_interp_accel_free (P2acc);
  gsl_interp_accel_free (P3acc);

  for (int i = 0; i <= Ns; i++)
    {
      gsl_spline_free (HHspline[i]);
      gsl_spline_free (VVspline[i]);
      gsl_spline_free (HPspline[i]);
      gsl_spline_free (VPspline[i]);

      gsl_interp_accel_free (HHacc[i]);
      gsl_interp_accel_free (VVacc[i]);
      gsl_interp_accel_free (HPacc[i]);
      gsl_interp_accel_free (VPacc[i]);
    }
  delete[] HHspline; delete[] VVspline; delete[] HPspline; delete[] VPspline;
  delete[] HHacc;    delete[] VVacc;    delete[] HPacc;    delete[] VPacc;

  delete[] cmu; delete[] ceta; delete[] seta; delete[] eeta; delete[] R2grgz; delete[] R2grge;  

  gsl_spline_free (Rrzspline);
  gsl_spline_free (Rrespline);
  gsl_spline_free (Rbspline);
  gsl_spline_free (Zbspline);

  gsl_interp_accel_free (Rrzacc);
  gsl_interp_accel_free (Rreacc);
  gsl_interp_accel_free (Rbacc);
  gsl_interp_accel_free (Zbacc);

  delete[] Rcoil; delete[] Zcoil;   delete[] Icoil;  delete[] Psix;
  delete[] Xi;    delete[] Upsilon; delete[] Lambda; delete[] Chi;

  delete[] tbound; delete[] Rbound; delete[] Zbound;
  delete[] dRdthe; delete[] dZdthe;

  delete[] mres; delete[] qres;  delete[] rres;   delete[] qerr;
  delete[] sres; delete[] DIres; delete[] nuLres; delete[] nuSres;
  delete[] Jres; delete[] DRres;

  delete[] S13res; delete[] taures; delete[] ieres; delete[] QEres;
  delete[] Qeres;  delete[] Qires;  delete[] Dres;  delete[] Pmres;
  delete[] Peres;  delete[] Teres;  delete[] Dcres;

  delete[] MPOL; delete[] mpol;
     
  delete[] Rgrid; delete[] hode; delete[] eode;

  delete[] Fval; delete[] Wval; delete[] deltaW;

  delete[] rf;
}  

// ##########################
// Function to read JSON file
// ##########################
json TJ::ReadJSONFile (const string& filename)
{
  ifstream JSONFile (filename);
  json     JSONData;

  if (JSONFile.is_open ())
    {
      try
	{
	  JSONFile >> JSONData;
        }
      catch (json::parse_error& e)
	{
	  cerr << "Unable to parse JSON file: " << e.what() << endl;
	  exit (1);
        }
      JSONFile.close ();
    }
  else
    {
      cerr << "Unable to open JSON file: " << filename << endl;
      exit (1);
    }

  return JSONData;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* TJ::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("TJ::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* TJ::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("TJ::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

