// TJ.cpp

#include "TJ.h"

// ###########
// Constructor
// ###########
TJ::TJ ()
{
  // -----------------------------------------
  // Ensure that directory ../Outputs/TJ exits
  // -----------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/TJ"))
    {
      exit (1);
    }

  // ------------------------------------------------
  // Ensure that directory ../Outputs/WriteEFIT exits
  // ------------------------------------------------
  if (!CreateDirectory ("../Outputs/WriteEFIT"))
    {
      exit (1);
    }

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

  // ---------------------------------------
  // Read control parameters from JSON files
  // ---------------------------------------
  string JSONFilename = "../Inputs/TJ.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR    = JSONData["NTOR"].get<int>();
  MMIN    = JSONData["MMIN"].get<int>();
  MMAX    = JSONData["MMAX"].get<int>();

  EQLB    = JSONData["EQLB"] .get<int>();
  FREE    = JSONData["FREE"] .get<int>();
  FVAL    = JSONData["FVAL"] .get<int>();
  RMP     = JSONData["RMP"]  .get<int>();
  IDEAL   = JSONData["IDEAL"].get<int>();
  XI      = JSONData["XI"]   .get<int>();
  INTR    = JSONData["INTR"] .get<int>();
  RWM     = JSONData["RWM"]  .get<int>();

  EPS     = JSONData["EPS"]  .get<double>();
  DEL     = JSONData["DEL"]  .get<double>();
  NFIX    = JSONData["NFIX"] .get<int>   ();
  NDIAG   = JSONData["NDIAG"].get<int>   ();

  NULC    = JSONData["NULC"]   .get<double>();
  ITERMAX = JSONData["ITERMAX"].get<int>   ();

  acc     = JSONData["acc"] .get<double>();
  h0      = JSONData["h0"]  .get<double>();
  hmin    = JSONData["hmin"].get<double>();
  hmax    = JSONData["hmax"].get<double>();

  EPSF    = JSONData["EPSF"].get<double>();

  JSONFilename = "../Inputs/Layer.json";
  JSONData     = ReadJSONFile (JSONFilename);
  
  B0    = JSONData["B0"]   .get<double>();
  R0    = JSONData["R0"]   .get<double>();
  n0    = JSONData["n0"]   .get<double>();
  alpha = JSONData["alpha"].get<double>();
  Zeff  = JSONData["Zeff"] .get<double>();
  Mion  = JSONData["Mion"] .get<double>();
  Chip  = JSONData["Chip"] .get<double>();
  Teped = JSONData["Teped"].get<double>();

  // ------------
  // Sanity check
  // ------------
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

  // No resonant magnetic perturbation calculation for fixed boundary
  if (FREE < 0)
    RMP = 0;

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\nClass TJ::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Calculation parameters:\n");
  printf ("ntor = %3d        mmin  = %3d        mmax = %3d        eps     = %10.3e del  = %10.3e\n",
	  NTOR, MMIN, MMAX, EPS, DEL);
  printf ("nfix = %3d        ndiag = %3d       nulc = %10.3e itermax = %3d        EQLB = %1d FREE = %1d FVAL = %1d RMP = %1d IDEAL = %1d XI = %1d INTR = %1d RWM = %1d\n",
	  NFIX, NDIAG, NULC, ITERMAX, EQLB, FREE, FVAL, RMP, IDEAL, XI, INTR, RWM);
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
  if (EQLB)
    return;
 
  // Calculate metric data at plasma boundary
  CalculateMetricBoundary ();

  // Calculate vacuum matrices
  GetVacuumBoundary ();

  // Calculate wall matrices
  GetVacuumWall ();

  // Find rational surfaces
  FindRational ();

  // Calculate resonant layer data
  GetLayerData ();

  // Solve outer region odes
  ODESolve ();

  // Determine tearing mode dispersion relation and tearing eigenfunctions
  FindDispersion ();

  // Calculate unreconnected eigenfunction visualization data
  VisualizeEigenfunctions ();

  if (RMP)
    {
      // Calculate resonant magnetic perturbation data
      CalculateResonantMagneticPerturbation ();

      // Calculate resonant magnetic perturbation respose visualization data
      VisualizeRMP ();
    }

  // Calculate ideal stability
  if (IDEAL)
    {
      CalculateIdealStability ();

      // Calculate resistive wall mode stability
      if (RWM)
	CalculateRWMStability ();
    }
  
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
  
  delete[] rr;   delete[] pp;  delete[] ppp; delete[] q; 
  delete[] s;    delete[] s2;  delete[] S1;  delete[] P1;
  delete[] P2;   delete[] P3;  delete[] g2;  delete[] p2;
  delete[] PsiN; delete[] S2;  delete[] S3;  delete[] s0;
  delete[] f;    delete[] Psi;

  gsl_spline_free (Pspline);
  gsl_spline_free (fspline);
  gsl_spline_free (g2spline);
  gsl_spline_free (p2spline);
  gsl_spline_free (ppspline);
  gsl_spline_free (pppspline);
  gsl_spline_free (qspline);
  gsl_spline_free (sspline);
  gsl_spline_free (s2spline);
  gsl_spline_free (s0spline);
  gsl_spline_free (S1spline);
  gsl_spline_free (S2spline);
  gsl_spline_free (S3spline);
  gsl_spline_free (P1spline);
  gsl_spline_free (P2spline);
  gsl_spline_free (P3spline);

  gsl_interp_accel_free (Pacc);
  gsl_interp_accel_free (facc);
  gsl_interp_accel_free (g2acc);
  gsl_interp_accel_free (p2acc);
  gsl_interp_accel_free (ppacc);
  gsl_interp_accel_free (pppacc);
  gsl_interp_accel_free (qacc);
  gsl_interp_accel_free (sacc);
  gsl_interp_accel_free (s2acc);
  gsl_interp_accel_free (s0acc);
  gsl_interp_accel_free (S1acc);
  gsl_interp_accel_free (S2acc);
  gsl_interp_accel_free (S3acc);
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

  delete[] cmu;  delete[] ceta;  delete[] seta;  delete[] eeta;  delete[] R2grgz;  delete[] R2grge;  

  gsl_spline_free (Rrzspline);
  gsl_spline_free (Rrespline);
  gsl_spline_free (Rbspline);
  gsl_spline_free (Zbspline);

  gsl_interp_accel_free (Rrzacc);
  gsl_interp_accel_free (Rreacc);
  gsl_interp_accel_free (Rbacc);
  gsl_interp_accel_free (Zbacc);

  delete[] Rcoil;   delete[] Zcoil;   delete[] Icoil;  delete[] Psix;
  delete[] Xi;

  if (RMP)
    {
      delete[] Upsilon; delete[] Lambda; delete[] Chi; delete[] Psirmps; delete[] Psixs;
    }
  
  delete[] tbound; delete[] Rbound; delete[] Zbound;
  delete[] dRdthe; delete[] dZdthe;

  delete[] wwall;  delete[] Rwall;  delete[] Zwall;
  
  delete[] mres; delete[] qres;  delete[] rres;   delete[] qerr;
  delete[] sres; delete[] DIres; delete[] nuLres; delete[] nuSres;
  delete[] Jres; delete[] DRres; delete[] Flarge; delete[] Fsmall;
  delete[] Pres;

  delete[] S13res; delete[] taures; delete[] ieres; delete[] QEres;
  delete[] Qeres;  delete[] Qires;  delete[] Dres;  delete[] Pmres;
  delete[] Peres;  delete[] Teres;  delete[] Dcres;

  delete[] MPOL; delete[] mpol;
     
  delete[] Rgrid; delete[] hode; delete[] eode; delete[] Pgrid;

  delete[] Fval;

  if (IDEAL)
    {
      delete[] Uval;    delete[] deltaW; delete[] deltaWv;
      delete[] deltaWp; delete[] gammax; delete[] gamma;
      delete[] Wval;    delete[] Vval;
    }

  if (RWM)
    {
      delete[] FFvl;
    }

  delete[] rf; delete[] rho;
}

// ########################################
// Function to strip comments from a string
// ########################################
string TJ::stripComments (const string& input)
{
  stringstream result;
  bool         inSingleLineComment = false;
  bool         inMultiLineComment  = false;

  for (size_t i = 0; i < input.size(); ++i)
    {
      // Start of single-line comment (//)
      if (!inMultiLineComment && input[i] == '/' && i + 1 < input.size() && input[i + 1] == '/')
	{
	  inSingleLineComment = true;
	  i++; 
	}
      // Start of multi-line comment (/* ... */)
      else if (!inSingleLineComment && input[i] == '/' && i + 1 < input.size() && input[i + 1] == '*')
	{
	  inMultiLineComment = true;
	  i++; 
	}
      // End of single-line comment
      else if (inSingleLineComment && input[i] == '\n')
	{
	  inSingleLineComment = false;
	  result << input[i];
	}
      // End of multi-line comment
      else if (inMultiLineComment && input[i] == '*' && i + 1 < input.size() && input[i + 1] == '/')
	{
	  inMultiLineComment = false;
	  i++; 
	}
      // Regular characters outside comments
      else if (!inSingleLineComment && !inMultiLineComment)
	{
	  result << input[i];
	}
    }
  
  return result.str();
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
	  // Strip any comments from JSON file
	  stringstream buffer;
	  buffer << JSONFile.rdbuf ();
	  JSONData = json::parse (stripComments (buffer.str ()));
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

// ################################################################
// Function to check that directory exists, and create it otherwise
// ################################################################
bool TJ::CreateDirectory (const char* path)
{
  struct stat st = {0};
  
  if (stat (path, &st) == -1)
    {
#ifdef _WIN32
      if (mkdir (path) != 0)
	{
	  printf ("Error creating directory: %s\n", path);
	  return false;
	}
#else
      if (mkdir (path, 0700) != 0)
	{
	  printf ("Error creating directory: %s\n", path);
	  return false;
	}
#endif
    }
  
  return true;
}



