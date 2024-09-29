// Plasma.cpp

// PROGRAM ORGANIATION:
//
// void Flux:: Stage2ReadData          ()
// void Flux:: Stage2CalcQ             ()
// void Flux:: Stage2CalcStraightAngle ()
// void Flux:: Stage2CalcBoundary      ()

#include "Flux.h"

// #############################################
// Function to read data for Stage2 calculations
// #############################################
void Flux::Stage2ReadData ()
{
  printf ("Reading data from EFIT file:\n");
  fflush (stdout);
  
  // ..............
  // Read R0 and B0
  // ..............
  FILE* file = OpenFiler ((char*) "Outputs/Flux/R0B0.txt");
  if (fscanf (file, "%lf %lf", &R0, &B0) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/R0B0.txt\n");
      exit (1);
    }
  fclose (file);
  
  // ......................
  // Read bounding box data
  // ......................
  file = OpenFiler ((char*) "Outputs/Flux/Box.txt");
  if (fscanf (file, "%lf %lf %lf %lf", &RLEFT, &ZLOW, &RRIGHT, &ZHIGH) != 4)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Box.txt\n");
      exit (1);
    }
  fclose (file);

  // .......................
  // Read magnetic axis data
  // .......................
  file = OpenFiler ((char*) "Outputs/Flux/Axis.txt");
  if (fscanf (file, "%lf %lf", &Raxis, &Zaxis) != 2)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Axis.txt\n");
      exit (1);
    }
  fclose (file);

  // ................
  // Read array sizes
  // ................
  file = OpenFiler ((char*) "Outputs/Flux/Points.txt");
  if (fscanf (file, "%d %d %d %d", &NRPTS, &NZPTS, &NBPTS, &NLPTS) != 4)
    {
      printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Points.txt\n");
      exit (1);
    }
  fclose (file);

  // ..............
  // Read R, Z grid
  // ..............
  RPTS = new double[NRPTS];  // R array
  ZPTS = new double[NZPTS];  // Z array

  file = OpenFiler ((char*) "Outputs/Flux/R.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf", &RPTS[i]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/R.txt\n");
	exit (1);
      }
  fclose (file);
  file = OpenFiler ((char*) "Outputs/Flux/Z.txt");
  for (int j = 0; j < NZPTS; j++)
    if (fscanf (file, "%lf", &ZPTS[j]) != 1)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Z.txt\n");
	exit (1);
      }
  fclose (file);

  // ..............................
  // Read boundary and limiter data
  // ..............................
  RBPTS = new double[NBPTS];  // R coordinates of boundary
  ZBPTS = new double[NBPTS];  // Z coordinates of boundary

  file = OpenFiler ((char*) "Outputs/Flux/Boundary.txt");
  for (int i = 0; i < NBPTS; i++)
    if (fscanf (file, "%lf %lf", &RBPTS[i], &ZBPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Boundary.txt\n");
	exit (1);
      }
  fclose (file);

  RLPTS = new double[NLPTS];  // R coordinates of limiter
  ZLPTS = new double[NLPTS];  // Z coordinates of limiter

  file = OpenFiler ((char*) "Outputs/Flux/Limiter.txt");
  for (int i = 0; i < NLPTS; i++)
    if (fscanf (file, "%lf %lf", &RLPTS[i], &ZLPTS[i]) != 2)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Limiter.txt\n");
	exit (1);
      }
  fclose (file);
  printf ("R0    = %10.3e  B0     = %10.3e\n", R0,  B0);
  printf ("RLEFT = %10.3e  RRIGHT = %10.3e  ZLOW = %10.3e  ZHIGH = %10.3e\n",
	   RLEFT, RRIGHT, ZLOW, ZHIGH);

  // ........
  // Read Psi
  // ........
  PSIARRAY.resize (NRPTS, NZPTS);  // Psi (R, Z)

  double val; int ival;
  file = OpenFiler ((char*) "Outputs/Flux/PsiSequential.txt");
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {	
	if (fscanf (file, "%d %d %lf", &ival, &ival, &val) != 3)
	  {
	    printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/PsiSequential.txt\n");
	    exit (1);
	  }
	PSIARRAY (i, j) = val;
      }
  fclose (file);

  // ....................................................
  // Modify PSI such that PsiBoundary = 0 and PsiAxis = 1
  // ....................................................
  double PsiAxis     = Interpolate (Raxis,    Zaxis,    NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 0, 0);
  double PsiBoundary = Interpolate (RBPTS[0], ZBPTS[0], NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 0, 0);
  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) - PsiBoundary;
	PSIARRAY (i, j) = val;
      }

  PsiAxis     = Interpolate (Raxis,    Zaxis,    NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 0, 0);
  PsiBoundary = Interpolate (RBPTS[0], ZBPTS[0], NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 0, 0);
  Psic        = PsiAxis / (R0*R0*B0);

  for (int i = 0; i < NRPTS; i++)
    for (int j = 0; j < NZPTS; j++)
      {
	double val = PSIARRAY (i, j) /PsiAxis;
	PSIARRAY (i, j) = val;
      }

  printf ("PsiAxis = %10.3e  PsiBoundary = %10.3e  PsiAxis /(R0*R0*B0) = %10.3e\n", PsiAxis, PsiBoundary, Psic);

  // .............................
  // Read equilibrium profile data
  // .............................
  PSIN = new double[NRPTS]; // PSI_N array 
  G    = new double[NRPTS]; // g
  Pr   = new double[NRPTS]; // p
  GGp  = new double[NRPTS]; // g dg/dpsi
  Prp  = new double[NRPTS]; // dp/dpsi
  Q    = new double[NRPTS]; // q

  file = OpenFiler ((char*) "Outputs/Flux/Profiles.txt");
  for (int i = 0; i < NRPTS; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf", &PSIN[i], &G[i], &Pr[i], &GGp[i], &Prp[i], &Q[i]) != 6)
      {
	printf ("FLUX::Stage2ReadData: Error reading Outputs/Flux/Profiles.txt\n");
	exit (1);
      }
  fclose (file);
}

// ######################################
// Function to calculate Stage2 q profile
// ######################################
void Flux::Stage2CalcQ ()
{
  // ....................................................
  // Find closest equilibrium grid-point to magnetic axis
  // ....................................................
  double rmin = 1.e6;
  for (int i = 0; i < NRPTS; i++)
    if (fabs (RPTS[i] - Raxis) < rmin)
      {
	rmin = fabs (RPTS[i] - Raxis);
	ic   = i;
      }
  rmin = 1.e6;
  for (int j = 0; j < NZPTS; j++)
    if (fabs (ZPTS[j] - Zaxis) < rmin)
      {
	rmin = fabs (ZPTS[j] - Zaxis);
	jc   = j;
      }

  // ..........................................................................
  // Find closest equilibrium grid-point to inner magnetic boundary on midplane
  // ..........................................................................
  ia = 0;
  for (int i = 0; i <= ic; i++)
    {
      if (PSIARRAY (i, jc) < 0.)
	ia++;
    }
  Rbound = (RPTS[ia-1]*PSIARRAY (ia, jc) - RPTS[ia]*PSIARRAY (ia-1, jc))
    /(PSIARRAY (ia, jc) - PSIARRAY (ia-1, jc));
 
  L = ic - ia + 2; // Number of points in Psi (R, Zaxis) array

  // .................................................
  // Calculate PSI array at magnetic axis level (in Z)
  // .................................................
  PSIAXIS = new double[NZPTS];

  double dZ = ZPTS[1] - ZPTS[0];
  for (int i = 0; i < NZPTS; i++)
    {
      if (ZPTS[jc] > Zaxis)
	PSIAXIS[i] = (PSIARRAY (i, jc-1) * (ZPTS[jc] - Zaxis)   + PSIARRAY (i, jc)   * (Zaxis - ZPTS[jc-1])) /dZ;
      else
	PSIAXIS[i] = (PSIARRAY (i, jc)   * (ZPTS[jc+1] - Zaxis) + PSIARRAY (i, jc+1) * (Zaxis - ZPTS[jc])) /dZ;
    }
  
  // ............................
  // Setup Psi (R, Zaxis) profile
  // ............................
  s = new double[L]; // Array of s = sqrt (1 - PSIAXIS) values

  for (int l = 0; l < L-2; l++)
    s[L-2-l] = sqrt (1. - PSIAXIS[l+ia]);
  s[0]   = 0.;
  s[L-1] = 1.;

  // ...........................
  // Setup R (Z = Zaxis) profile
  // ...........................
  Rs = new double[L]; // Array of R(s) values

  for (int l = 0; l < L-2; l++)
    Rs[L-2-l] = RPTS[l+ia];
  Rs[0]   = Raxis;
  Rs[L-1] = Rbound;

  // ......................
  // Set up Stage2 Psi grid
  // ......................
  P    = new double[NPSI];  // 1 - Psi array
  RP   = new double[NPSI];  // R(Psi) on inboard midplane
  rP   = new double[NPSI];  // r(Psi)
  GP   = new double[NPSI];  // g(Psi)
  QGP  = new double[NPSI];  // q(Psi)/g(Psi)
  QP   = new double[NPSI];  // q(Psi)
  FP   = new double[NPSI];  // f(Psi)
  PP   = new double[NPSI];  // P(Psi)
  GPP  = new double[NPSI];  // dg/dPsi
  PPP  = new double[NPSI];  // dP/dPsi
  GPX  = new double[NPSI];  // g dg/dPsi (check)
  PPX  = new double[NPSI];  // dP/dPsi (check)
  S    = new double[NPSI];  // sqrt (1 - Psi)
  QX   = new double[NPSI];  // q(Psi) from gFile
  PsiN = new double[NPSI];  // PsiN array
  
  for (int j = 0; j < NPSI; j++)
    {
      double ss = double (j+1) /double (NPSI);

      PsiN[j] = (1. - pow (1. - ss, PACK )) * (1. - pow (1. - ss, PACK ));
      P   [j] = 1. - PsiN[j];
      S   [j] = sqrt (1. - P[j]);
    }

  // ...............................
  // Calculate Stage2 R(Psi) profile
  // ...............................
  for (int j = 0; j < NPSI; j++)
    RP[j] = Interpolate (L, s, Rs, S[j], 0);

  // ....................................
  // Recalculate Stage2 P and PsiN arrays
  // ....................................
  for (int j = 0; j < NPSI; j++)
    {
      P[j]    = GetPsi (RP[j], Zaxis);
      PsiN[j] = 1. - P[j];
    }

  // .......................................
  // Calculate Stage2 g(Psi), P(Psi) profile
  // .......................................
  for (int j = 0; j < NPSI; j++)
    {
      GP [j] = Interpolate (NRPTS, PSIN, G,    PsiN[j], 0);
      PP [j] = Interpolate (NRPTS, PSIN, Pr,   PsiN[j], 0);
      GPP[j] = Interpolate (NRPTS, PSIN, GGp,  PsiN[j], 0) /GP[j];
      PPP[j] = Interpolate (NRPTS, PSIN, Prp,  PsiN[j], 0);
      QX [j] = Interpolate (NRPTS, PSIN, Q,    PsiN[j], 0);
     }

  for (int j = 0; j < NPSI; j++)
    {
      GPX[j] = - GP[j] * Interpolate (NRPTS, PSIN, G,  PsiN[j], 1) /Psic;
      PPX[j] = -         Interpolate (NRPTS, PSIN, Pr, PsiN[j], 1) /Psic;
    }

  // ......................................
  // Calculate Stage2 q(Psi)/g(Psi) profile
  // ......................................
  printf ("Calculating q(Psi)/g(Psi) profile:\n");
  fflush (stdout);
  CalcQGP ();

  // .............................
  // Calculate Stage2 r(P) profile
  // .............................
  printf ("Calculating r(Psi) profile:\n");
  fflush (stdout);
  CalcrP ();

  for (int j = 0; j < NPSI; j++)
    {
      FP[j] = rP[j] * GP[j] /QP[j];
    }

  // ...........
  // Find qa, ra
  // ...........
  qa = QP[NPSI-1];
  ra = rP[NPSI-1];
  printf ("qa = %10.3e  ra = %10.3e  a = %10.3e (m)\n",
	  qa, ra, ra*R0);
}

// #################################################
// Calculate Stage2 straight angle coordinate system
// #################################################
void Flux::Stage2CalcStraightAngle ()
{
  // ..................
  // Set up theta array
  // ..................
  th = new double[NTHETA]; 
  for (int k = 0; k < NTHETA; k++)
    {
      double t = double (k) /double (NTHETA-1);
      th[k]    = 2.* M_PI * t;
    }

  // ...............
  // Allocate memory
  // ...............
  RRst = gsl_matrix_alloc (NPSI, NTHETA);
  ZZst = gsl_matrix_alloc (NPSI, NTHETA);
  RRr  = gsl_matrix_alloc (NPSI, NTHETA);
  RRth = gsl_matrix_alloc (NPSI, NTHETA);
  ZZr  = gsl_matrix_alloc (NPSI, NTHETA);
  ZZth = gsl_matrix_alloc (NPSI, NTHETA);
  Jac  = gsl_matrix_alloc (NPSI, NTHETA);
  Jax  = gsl_matrix_alloc (NPSI, NTHETA);

  // ..........................................
  // ..........................................
  printf ("Calculating straight angle coordinate system:\n");
  fflush (stdout);
  CalcStraightAngleGeneral ();

  // .............................
  // Calculate partial derivatives
  // .............................
  for (int j = 0; j < NPSI; j++)
    {
      qgp = fabs (Psic) * QP[j] /GP[j];

      for (int k = 0; k < NTHETA; k++)
	{
	  double Rval = gsl_matrix_get (RRst, j, k);
	  double Zval = gsl_matrix_get (ZZst, j, k);
	  double PsiR = GetPsiR (Rval, Zval);
	  double PsiZ = GetPsiZ (Rval, Zval);

	  gsl_matrix_set (RRth, j, k, - qgp * Rval * PsiZ /rP[j]);
	  gsl_matrix_set (ZZth, j, k,   qgp * Rval * PsiR /rP[j]);
	}
    }

  for (int k = 0; k < NTHETA; k++)
    {
      for (int j = 0; j < NPSI; j++)
	{
	  double dRdr = Interpolate1 (NPSI, rP, RRst, k, rP[j], 1);
	  double dZdr = Interpolate1 (NPSI, rP, ZZst, k, rP[j], 1);

	  gsl_matrix_set (RRr, j, k, dRdr);
	  gsl_matrix_set (ZZr, j, k, dZdr);
	}
    }

  // Calculate Jacobian (numerical and analytic)
  for (int j = 0; j < NPSI; j++)
    for (int k = 0; k < NTHETA; k++)
      {
	double Rval = gsl_matrix_get (RRst, j, k);
	double Jval = gsl_matrix_get (RRth, j, k) * gsl_matrix_get (ZZr, j, k)
	  - gsl_matrix_get (ZZth, j, k) * gsl_matrix_get (RRr, j, k);

	gsl_matrix_set (Jac, j, k, Jval);
	gsl_matrix_set (Jax, j, k, Rval);
      }
}

// ##############################
// Calculate Stage2 boundary data
// ##############################
void Flux::Stage2CalcBoundary ()
{
  Rbndry   = new double[NTHETA];
  Zbndry   = new double[NTHETA];
  dRdtheta = new double[NTHETA];
  dZdtheta = new double[NTHETA];

  for (int k = 0; k < NTHETA; k++)
    {
      Rbndry  [k] = gsl_matrix_get (RRst, NPSI-1, k);
      Zbndry  [k] = gsl_matrix_get (ZZst, NPSI-1, k);
      dRdtheta[k] = gsl_matrix_get (RRth, NPSI-1, k) * ra;
      dZdtheta[k] = gsl_matrix_get (ZZth, NPSI-1, k) * ra;
    }
}
