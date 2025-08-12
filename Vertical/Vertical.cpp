// Vertical.cpp

#include "Vertical.h"

// ###########
// Constructor
// ###########
Vertical::Vertical ()
{
  // -----------------------------------------------
  // Ensure that directory ../Outputs/Vertical exits
  // -----------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Vertical"))
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
  string JSONFilename = "../Inputs/Vertical.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  MMIN    = JSONData["MMIN"].get<int>();
  MMAX    = JSONData["MMAX"].get<int>();

  EQLB    = JSONData["EQLB"].get<int>();
  
  EPS     = JSONData["EPS"]  .get<double>();
  NFIX    = JSONData["NFIX"] .get<int>   ();
  NDIAG   = JSONData["NDIAG"].get<int>   ();

  acc     = JSONData["acc"] .get<double>();
  h0      = JSONData["h0"]  .get<double>();
  hmin    = JSONData["hmin"].get<double>();
  hmax    = JSONData["hmax"].get<double>();

  EPSF    = JSONData["EPSF"].get<double>();

  JSONFilename = "../Inputs/Equilibrium.json";
  JSONData     = ReadJSONFile (JSONFilename);
  
  SRC   = JSONData["SRC"]  .get<int>   ();
  B0    = JSONData["B0"]   .get<double>();
  R0    = JSONData["R0"]   .get<double>();
  n0    = JSONData["n0"]   .get<double>();
  alpha = JSONData["alpha"].get<double>();
  Zeff  = JSONData["Zeff"] .get<double>();
  Mion  = JSONData["Mion"] .get<double>();
  Chip  = JSONData["Chip"] .get<double>();
  Teped = JSONData["Teped"].get<double>();
  neped = JSONData["neped"].get<double>();

  JSONFilename = "../Inputs/TJ.json";
  JSONData     = ReadJSONFile (JSONFilename);
  VIZ          = JSONData["VIZ"].get<int> ();

  // ------------
  // Sanity check
  // ------------
  if (MMAX < MMIN)
    {
      printf ("Vertical: Error - MMIN must be less that MMAX\n");
      exit (1);
    }
  if (MMAX < MMIN)
    {
      printf ("Vertical: Error - MMIN must be less that MMAX\n");
      exit (1);
    }
  if (EPS <= 0.)
    {
      printf ("Vertical: Error - EPS must be positive\n");
      exit (1);
    }
  if (NFIX < 0)
    {
      printf ("Vertical: Error - NFIX cannot be negative\n");
      exit (1);
    }
  if (NDIAG < 0)
    {
      printf ("Vertical: Error - NDIAG cannot be less that two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("Vertical:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("Vertical:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("Vertical:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("Vertical:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("Vertical:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (EPSF <= 0.)
    {
      printf ("Vertical: Error - EPSF must be positive\n");
      exit (1);
    }
  if (R0 <= 0.)
    {
      printf ("Vertical: Error - R0 must be positive\n");
      exit (1);
    }
  if (n0 <= 0.)
    {
      printf ("Vertical: Error - n0 must be positive\n");
      exit (1);
    }
  if (Zeff < 1.)
    {
      printf ("Vertical: Error - Zeff cannot be less than unity\n");
      exit (1);
    }
  if (Mion < 1.)
    {
      printf ("Vertical: Error - Mion cannot be less than unity\n");
      exit (1);
    }
  if (Chip <= 0.)
    {
      printf ("Vertical: Error - Chip must be positive\n");
      exit (1);
    }
  if (Teped <= 0.)
    {
      printf ("Vertical: Error - Teped must be positive\n");
      exit (1);
    }
  if (neped <= 0.)
    {
      printf ("Vertical: Error - neped must be positive\n");
      exit (1);
    }
}

// ##########
// Destructor
// ##########
Vertical::~Vertical ()
{
}

// #########################
// Function to solve problem
// #########################
void Vertical::Solve ()
{
  // .............................................................................
  // Call class Equilibrium to construct aspect-ratio expanded tokamak equilibrium
  // .............................................................................
  {
    Equilibrium equilibrium;

    if (!SRC)
      equilibrium.Setnu ();
    
    equilibrium.Solve ();
  }

  if (!EQLB)
    {
      // -----------------------------
      // Output calculation parameters
      // -----------------------------
      printf ("\nClass Vertical::\n");
      printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
      printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
      printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
      printf ("Calculation flags:\n");
      printf ("EQLB = %1d  VIZ = %1d\n",
	      EQLB, VIZ);
      printf ("Calculation parameters:\n");
      printf ("mmin = %3d        mmax  = %-3d        eps   = %10.3e\n",
	      MMIN, MMAX, EPS);
      printf ("nfix = %3d        ndiag = %3d\n",
	      NFIX, NDIAG);
      printf ("acc  = %10.3e h0    = %10.3e hmin  = %10.3e hmax    = %10.3e epsf = %10.3e\n",
	      acc, h0, hmin, hmax, EPSF);
      printf ("B0   = %10.3e R0    = %10.3e n0    = %10.3e alpha   = %10.3e Zeff = %10.3e\n",
	      B0, R0, n0, alpha, Zeff);
      printf ("Mion = %10.3e Chip  = %10.3e Teped = %10.3e neped   = %10.3e\n",
	      Mion, Chip, Teped, neped);

      // Set toroidal and poloidal mode numbers
      SetModeNumbers ();

      // Read equilibrium data
      ReadEquilibrium ();

      // Calculate metric data at plasma boundary
      CalculateMetricBoundary ();

      // Calculate vacuum matrices
      GetVacuumBoundary ();

      // Calculate wall matrices
      GetVacuumWall ();
        
      // Solve outer region odes
      ODESolve ();

      // Calculate ideal stability
      CalculateIdealStability ();
      
      if (VIZ)
	{
	  // Calculate ideal eigenfunction visualization data
	  VisualizeEigenfunctions ();
	}
    
      // Write program data to Netcdf file
      WriteNetcdf ();
      
      // Clean up
      CleanUp ();
    }
}

// ##################################################
// Function to set toroidal and poloidal mode numbers
// ##################################################
void Vertical::SetModeNumbers ()
{
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
void Vertical::CleanUp ()
{
  printf ("Cleaning up:\n");
  
  delete[] rr;   delete[] pp;  delete[] ppp; delete[] q; 
  delete[] s;    delete[] s2;  delete[] S1;  delete[] P1;
  delete[] P2;   delete[] P3;  delete[] g2;  delete[] p2;
  delete[] PsiN; delete[] S2;  delete[] S3;  delete[] s0;
  delete[] f;    delete[] Psi; delete[] nep; delete[] Tep;
  delete[] ne;   delete[] Te;  delete[] S4;  delete[] S5;
  delete[] Sig;

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
  gsl_spline_free (S4spline);
  gsl_spline_free (S5spline);
  gsl_spline_free (Sigspline);
  gsl_spline_free (P1spline);
  gsl_spline_free (P2spline);
  gsl_spline_free (P3spline);
  gsl_spline_free (nespline);
  gsl_spline_free (Tespline);
  gsl_spline_free (nepspline);
  gsl_spline_free (Tepspline);

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
  gsl_interp_accel_free (S4acc);
  gsl_interp_accel_free (S5acc);
  gsl_interp_accel_free (Sigacc);
  gsl_interp_accel_free (P1acc);
  gsl_interp_accel_free (P2acc);
  gsl_interp_accel_free (P3acc);
  gsl_interp_accel_free (neacc);
  gsl_interp_accel_free (Teacc);
  gsl_interp_accel_free (nepacc);
  gsl_interp_accel_free (Tepacc);

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
  
  delete[] rho;
}

