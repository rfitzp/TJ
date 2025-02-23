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
  VIZ     = JSONData["VIZ"]  .get<int>();
  IDEAL   = JSONData["IDEAL"].get<int>();
  XI      = JSONData["XI"]   .get<int>();
  INTR    = JSONData["INTR"] .get<int>();
  RWM     = JSONData["RWM"]  .get<int>();
  LAYER   = JSONData["LAYER"].get<int>();

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
  printf ("Calculation flags:\n");
  printf ("EQLB = %1d FREE = %1d FVAL = %1d RMP = %1d VIZ = %1d IDEAL = %1d XI = %1d INTR = %1d RWM = %1d LAYER = %1d\n",
	  EQLB, FREE, FVAL, RMP, VIZ, IDEAL, XI, INTR, RWM, LAYER);
  printf ("Calculation parameters:\n");
  printf ("ntor = %3d        mmin  = %3d        mmax = %3d        eps     = %10.3e del  = %10.3e\n",
	  NTOR, MMIN, MMAX, EPS, DEL);
  printf ("nfix = %3d        ndiag = %3d       nulc = %10.3e itermax = %3d\n",
	  NFIX, NDIAG, NULC, ITERMAX);
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

  if (VIZ)
    {
      // Calculate unreconnected eigenfunction visualization data
      VisualizeEigenfunctions ();
    }

  if (RMP)
    {
      // Calculate resonant magnetic perturbation data
      CalculateResonantMagneticPerturbation ();

      if (VIZ)
	{
	  // Calculate resonant magnetic perturbation respose visualization data
	  VisualizeRMP ();
	}
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
      delete[] FFvl; delete[] fw;
    }

  if (VIZ)
    delete[] rf;
  delete[] rho;
}

