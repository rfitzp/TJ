// NonLinear.cpp

#include "NonLinear.h"

// ###########
// Constructor
// ###########
NonLinear::NonLinear (int _machine)
{
  // --------------------------------------
  // Set default values of class parameters
  // --------------------------------------

  machine = _machine;

  // Machine setup
  if (machine == 1)
    {
      // %%%
      // RFP
      // %%%
      
      // Equilibrium parameters
      epsa   = 0.25;
      epsb   = 1.01 * epsa;
      q0     = 0.14;
      beta0  = 6.e-2;
      alphas = 4.;
      nus    = 1.6;
      alphap = 4.;
      nup    = 2.;
      
      // Mode parameters
      mpol1 = 1.;
      ntor1 = 8.;
      mpol2 = 0.;
      ntor2 = 1.;
      mpol3 = mpol1 + mpol2;
      ntor3 = ntor1 + ntor2;
      W1    = 1.e-2;
      W2    = 1.e-2;
      W3    = 1.e-2;
    }
  else if (machine == 2)
    {
      // %%%%%%%
      // Tokamak
      // %%%%%%%
      
      // Equilibrium parameters
      epsa   = 0.25;
      epsb   = 1.1 * epsa;
      q0     = 0.9;
      beta0  = 0.1;
      alphas = 2.;
      nus    = 2.;
      alphap = 2.;
      nup    = 2.;
      
      // Mode parameters
      mpol1 = 1.;
      ntor1 = 1.;
      mpol2 = 2.;
      ntor2 = 1.;
      mpol3 = mpol1 + mpol2;
      ntor3 = ntor1 + ntor2;
      W1    = 1.e-3;
      W2    = 1.e-3;
      W3    = 1.e-3;
    }
  else
    exit(0);
  
  // Calculation parameters
  Ngrid  = 5000;
  Ncntr  = 500;
  eps    = 1.e-4;
  delta  = 1.e-7;
  TorCur = 0;

  // Adaptive integration parameters 
  acc     = 1.e-14;
  h0      = 1.e-6;
  hmin    = 1.e-10;
  hmax    = 1.e-4;
  maxrept = 50;
  flag    = 2;

  // ----------------------
  // Set RK4/RK5 parameters
  // ----------------------
  aa1 = 0.;
  aa2 = 1./5.;
  aa3 = 3./10.;
  aa4 = 3./5.;
  aa5 = 1.;
  aa6 = 7./8.;

  cc1 =  37./378.;
  cc3 = 250./621.;
  cc4 = 125./594.;
  cc6 = 512./1771.;

  ca1 = cc1 -  2825./27648.;
  ca3 = cc3 - 18575./48384.;
  ca4 = cc4 - 13525./55296.;
  ca5 =     -   277./14336.;
  ca6 = cc6 -     1./4.;

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

  // ---------------
  // Allocate memory
  // ---------------
  rr       = new double[Ngrid+1];
  ssigma   = new double[Ngrid+1];
  PP       = new double[Ngrid+1];
  PQ       = new double[Ngrid+1];
  BBphi    = new double[Ngrid+1];
  BBtheta  = new double[Ngrid+1];
  qq       = new double[Ngrid+1];
  qp       = new double[Ngrid+1];

  rrc      = new double[Ncntr+1];
  PPc      = new double[Ncntr+1];
  Btheta2j = new double[Ncntr+1];
  Bthetam  = new double[Ncntr+1];
  Bthetap  = new double[Ncntr+1];

  sssigma  = new double[Ngrid+1];
  PPP      = new double[Ngrid+1];
  BBBphi   = new double[Ngrid+1];
  BBBtheta = new double[Ngrid+1];
  qqq      = new double[Ngrid+1];

  Psi1     = new double[Ngrid+1];
  Psip1    = new double[Ngrid+1];
  Psi2     = new double[Ngrid+1];
  Psip2    = new double[Ngrid+1];
  Psi3     = new double[Ngrid+1];
  Psip3    = new double[Ngrid+1];

  JJ1      = new double[Ncntr+1];
  JJ2      = new double[Ncntr+1];
  JJ3      = new double[Ncntr+1];

  PPsi1    = new double[Ngrid+1];
  PPsip1   = new double[Ngrid+1];
  PPsi2    = new double[Ngrid+1];
  PPsip2   = new double[Ngrid+1];
  PPsi3    = new double[Ngrid+1];
  PPsip3   = new double[Ngrid+1];

  Xi1      = new double[Ngrid+1];
  Xi2      = new double[Ngrid+1];
  Xi3      = new double[Ngrid+1];

  Tau      = new double[Ngrid+1];
  Taup     = new double[Ngrid+1];

  Over     = new double[Ngrid+1];
  Overp    = new double[Ngrid+1];

  Atau     = gsl_interp_accel_alloc ();
  Ataup    = gsl_interp_accel_alloc ();
  APQ      = gsl_interp_accel_alloc ();

  Stau     = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);
  Staup    = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);
  SPQ      = gsl_spline_alloc (gsl_interp_cspline, Ngrid+1);

  // ------------------
  // Set up radial grid
  // ------------------
  double dr = 1. /double (Ngrid);
  for (int i = 0; i <= Ngrid; i++)
    rr[i] = double (i) * dr;

  // -----------------------
  // Set up control surfaces
  // -----------------------
  double drc = 1. /double (Ncntr);
  for (int i = 0; i <= Ncntr; i++)
    rrc[i] = double (i) * drc;
  rrc[Ncntr] += 1.e-12;
}

// ###########
// Destructor
// ###########
NonLinear::~NonLinear ()
{
  delete[] rr;      delete[] ssigma; delete[] PP;       delete[] BBphi;    delete[] BBtheta; delete[] qq;     delete[] qp;
  delete[] rrc;     delete[] PPc;    delete[] Btheta2j; delete[] Bthetam;  delete[] Bthetap;
  delete[] sssigma; delete[] PPP;    delete[] BBBphi;   delete[] BBBtheta; delete[] qqq;
  delete[] Psi1;    delete[] Psip1;  delete[] Psi2;     delete[] Psip2;    delete[] Psi3;    delete[] Psip3;
  delete[] PPsi1;   delete[] PPsip1; delete[] PPsi2;    delete[] PPsip2;   delete[] PPsi3;   delete[] PPsip3;
  delete[] JJ1;     delete[] JJ2;    delete[] JJ3;      delete[] Xi1;      delete[] Xi2;     delete[] Xi3;
  delete[] Tau;     delete[] Taup;   delete[] Over;     delete[] Overp;    delete[] PQ;

  gsl_spline_free       (Stau); gsl_spline_free       (Staup); gsl_spline_free       (SPQ);
  gsl_interp_accel_free (Atau); gsl_interp_accel_free (Ataup); gsl_interp_accel_free (APQ);
}

// ####################################
// Function to evaluate stability index
// ####################################
void NonLinear::GetDelta (int _mpol, int _ntor, double& _Delta, double& _Psirv)
{
  // Calculate continuous zero pressure equilibrium
  CalcZero ();
  
  // Calculate continuous pressure equilibrium
  CalcEquilibrium ();

  // Calculate stepped pressure equilibrium
  CalcStepped ();

  // Find rational surface
  double _mpol_ = double (_mpol);
  double _ntor_ = double (_ntor);
  double qs     = _mpol_ /_ntor_;
  double _rs    = Findrs (qs);
  if (_rs > 0.)
    {
    }
  else
    {
      printf ("\nmpol = %11.4e  ntor = %11.4e - No rational surface in plasma!\n", _mpol_, _ntor_);
      exit(1);
    }

  // Get control surfaces
  GetJ (_mpol_, _ntor_, _rs, W1, JJ1);
  
  // Calculate tearing eigenfuction 
  GetDelta (_mpol_, _ntor_, _rs, JJ1, _Delta, _Psirv);
}

// #########################
// Function to solve problem
// #########################
void NonLinear::Solve (int verbose, double& tau, double& taup, double& _Delta1, double& _Delta2, double& _Delta3)
{
  // Calculate continuous zero pressure equilibrium
  CalcZero ();
  
  // Calculate continuous pressure equilibrium
  CalcEquilibrium ();

  if (verbose)
    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nCalculating Equilibrium:\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  
  if (verbose)
    printf ("Zero pressure:       Theta = %11.4e  F = %11.4e\n", Theta, Frev);

  // Calculate stepped pressure equilibrium
  CalcStepped ();

  if (verbose)
    printf ("Continuous pressure: Theta = %11.4e  F = %11.4e\n", TTheta, FFrev);


  if (verbose)
    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nCalculating Continuous Eigenfunctions:\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

  // Find rational surface of first mode
  double qs = mpol1 /ntor1;
  rs1       = Findrs (qs);
  if (rs1 > 0.)
    {
      if (verbose)
	printf ("First mode:  mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol1, ntor1, rs1);
    }
  else
    {
      printf ("\nFirst mode:  mpol = %11.4e  ntor = %11.4e - No rational surface in plasma!\n", mpol1, ntor1);
      exit(1);
    }

  // Calculate tearing eigenfuction of first mode
  GetPsi (mpol1, ntor1, rs1, Delta1, Psi1, Psip1, verbose); 
  
  // Find rational surface of second mode
  qs = mpol2 /ntor2;
  rs2 = Findrs (qs);
  if (rs2 > 0.)
    {
      if (verbose)
	printf ("Second mode: mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol2, ntor2, rs2);
    }
  else
    {
      printf ("Second mode: mpol = %11.4e  ntor = %11.4e - No rational surface in plasma!\n", mpol2, ntor2);
      exit(1);
    }

  // Calculate tearing eigenfuction of second mode
  GetPsi(mpol2, ntor2, rs2, Delta2, Psi2, Psip2, verbose);

  // Find rational surface of third mode
  qs = mpol3 /ntor3;
  rs3 = Findrs (qs);
  if (rs3 > 0.)
    {
      if (verbose)
	printf ("Third mode:  mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol3, ntor3, rs3);
    }
  else
    {
      printf ("Third mode:  mpol = %11.4e  ntor = %11.4e - No rational surface in plasma!\n", mpol3, ntor3);
      exit(1);
    }
  
  // Calculate tearing eigenfuction of third mode
  GetPsi (mpol3, ntor3, rs3, Delta3, Psi3, Psip3, verbose);

  // Output rational surfaces
  FILE* file = OpenFilew ("Output/Rational.out");
  fprintf (file, "%e %e %e\n", rs1, rs2, rs3);
  fclose (file);

  // Output mode numbers
  file = OpenFilew ("Output/Mode.out");
  fprintf (file, "%d %d %d %d %d %d\n", int(mpol1), int(ntor1), int(mpol2), int(ntor2), int (mpol3), int(ntor3));
  fclose (file);

  // Output control surfaces
  file = OpenFilew ("Output/Control.out");
  for (int i = 1; i < Ncntr; i++)
    fprintf (file, "%e\n", rrc[i]);
  fclose (file);

  // Output tearing eigenfuctions
  file = OpenFilew ("Output/Eigenfunction.out");
  for (int i = 0; i <= Ngrid; i++)
    fprintf (file, "%e %e %e %e %e %e %e\n", rr[i], Psi1[i], Psip1[i], Psi2[i], Psip2[i], Psi3[i], Psip3[i]);
  fclose (file);

  // Calculate jump conditions at control surfaces for first, second, and third modes
  GetJ (mpol1, ntor1, rs1, W1, JJ1);
  GetJ (mpol2, ntor2, rs2, W2, JJ2);
  GetJ (mpol3, ntor3, rs3, W3, JJ3);

  if (verbose)
    printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nCalculating Stepped Eigenfunctions:\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  
  // Calculate stepped tearing eigenfunction of first mode
  if (verbose)
    printf ("First mode:  mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol1, ntor1, rs1);
  GetPsiStepped (mpol1, ntor1, rs1, W1, JJ1, DDelta1, PPsi1, PPsip1, Xi1, verbose);

  // Calculate stepped tearing eigenfunction of second mode
  if (verbose)
    printf ("Second mode: mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol2, ntor2, rs2);
  GetPsiStepped (mpol2, ntor2, rs2, W2, JJ2, DDelta2, PPsi2, PPsip2, Xi2, verbose);

  // Calculate stepped tearing eigenfunction of third mode
  if (verbose)
    printf ("Third mode:  mpol = %11.4e  ntor = %11.4e  rs = %11.4e\n", mpol3, ntor3, rs3);
  GetPsiStepped (mpol3, ntor3, rs3, W3, JJ3, DDelta3, PPsi3, PPsip3, Xi3, verbose);

  // Output stepped tearing eigenfuctions
  file = OpenFilew ("Output/EigenStepped.out");
  for (int i = 0; i <= Ngrid; i++)
    {
      fprintf (file, "%e %e %e %e %e %e %e %e %e %e\n",
	       rr[i], PPsi1[i], PPsip1[i], PPsi2[i], PPsip2[i], PPsi3[i], PPsip3[i], Xi1[i], Xi2[i], Xi3[i]);
    }
  fclose (file);

  _Delta1 = DDelta1;
  _Delta2 = DDelta2;
  _Delta3 = DDelta3;

  if (verbose)
    printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nCalculating Overlap Integrals:\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

  // Calculate overlap functions
  GetTau ();

  // Output overlap functions
  file = OpenFilew ("Output/Overlap.out");
  for (int i = 0; i <= Ngrid; i++)
    fprintf (file, "%e %e %e\n", rr[i], Tau[i], Taup[i]);
  fclose (file);

  // Set up interpolation of overlap functions
  gsl_spline_init (Stau,  rr, Tau,  Ngrid+1);
  gsl_spline_init (Staup, rr, Taup, Ngrid+1);

  // Calculate overlap integrals
  CalcOverlap ();
  tau  = Over[Ngrid];
  taup = Overp[Ngrid];

  // Output overlap integrals
  file = OpenFilew ("Output/Integral.out");
  for (int i = 0; i <= Ngrid; i++)
    fprintf (file, "%e %e %e\n", rr[i], Over[i], Overp[i]);
  fclose (file);
}

// ############################
// Function to set mode numbers
// ############################
void NonLinear::SetMode (int _n, int _k)
{
  mpol1 = 1.;
  ntor1 = double (_n);
  mpol2 = 0.;
  ntor2 = double (_k);
  mpol3 = mpol1 + mpol2;
  ntor3 = ntor1 + ntor2;
}

// #####################
// Function to set beta0
// #####################
void NonLinear::Setbeta0 (double _beta0)
{
  beta0 = _beta0;
}

// #################
// Function to set W
// #################
void NonLinear::SetW (double _W)
{
  W1 = _W;
  W2 = _W;
  W3 = _W;
}

// ##########################################################
// Function to calculate continuous zero pressure equilibrium
// ##########################################################
void NonLinear::CalcZero ()
{
  // Calculate parallel current and pressure profiles
  for (int i = 0; i <= Ngrid; i++)
    {
      ssigma[i] = GetSigma (rr[i]);
      PP    [i] = GetP     (rr[i]);
    }

  // Calculate magnetic field and safety-factor profiles
  BBphi[0]   = 1.;
  BBtheta[0] = 0.;
  qq[0]      = q0;
  PQ[0]      = 0.;

  double  r, h, t_err;
  int     rept; rhs_chooser = 0;
  double* y   = new double[4];
  double* err = new double[4];

  r     = eps;
  h     = h0;
  count = 0;

  y[0] = 1.;
  y[1] = (epsa /q0) * r;
  y[2] = 0.;
  y[3] = 0.;

  for (int i = 1; i <= Ngrid; i++)
    {
      do
	{
	  RK4RK5Adaptive (4, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      RK4RK5Fixed (4, r, y, err, rr[i] - r);

      BBphi[i]   = y[0];
      BBtheta[i] = y[1];
      qq[i]      = epsa * rr[i] * BBphi[i] /BBtheta[i];
      PQ[i]      = y[2];
    }

  Theta = y[1] /2./y[3];
  Frev  = y[0] /2./y[3];
  
  delete[] y; delete[] err;

  double PQa = PQ[Ngrid];
  for (int i = 0; i <= Ngrid; i++)
    PQ[i] -= PQa;

  // Set up interpolation of modified pressure profile
  gsl_spline_init (SPQ, rr, PQ, Ngrid+1);

  // Output equilibrium profiles
  FILE* file = OpenFilew ("Output/Zero.out");

  for (int i = 0; i <= Ngrid; i++)
    if (machine == 1)
      fprintf (file, "%e %e %e %e %e %e %e\n", rr[i], ssigma[i], PP[i], BBphi[i], BBtheta[i], qq[i], PQ[i]);
    else
      fprintf (file, "%e %e %e %e %e %e %e\n", rr[i], ssigma[i], PP[i], BBphi[i], BBtheta[i], qq[i], -PQ[i]);
      
  fclose (file);
}

// #####################################################
// Function to calculate continuous pressure equilibrium
// #####################################################
void NonLinear::CalcEquilibrium ()
{
  // Calculate parallel current and pressure profiles
  for (int i = 0; i <= Ngrid; i++)
    {
      ssigma[i] = GetSigma (rr[i]);
      PP    [i] = GetPQ    (rr[i]);
    }

  // Calculate magnetic field and safety-factor profiles
  BBphi[0]   = 1.;
  BBtheta[0] = 0.;
  qq[0]      = q0;
  qp[0]      = 0.;
 
  double  r, h, t_err;
  int     rept; rhs_chooser = 1;
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
	  RK4RK5Adaptive (3, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      RK4RK5Fixed (3, r, y, err, rr[i] - r);

      BBphi[i]   = y[0];
      BBtheta[i] = y[1];
      qq[i]      = epsa * rr[i] * BBphi[i] /BBtheta[i];
      qp[i]      = (2. - GetSigma (rr[i]) * qq[i] /epsa) * qq[i] /rr[i] - epsa * rr[i] * GetSigma (rr[i]);
    }

  TTheta = y[1] /2./y[2];
  FFrev  = y[0] /2./y[2];
  
  delete[] y; delete[] err;

  // Output equilibrium profiles
  FILE* file = OpenFilew ("Output/Equilibrium.out");

  for (int i = 0; i <= Ngrid; i++)
    if (machine == 1)
      fprintf (file, "%e %e %e %e %e %e %e\n", rr[i], ssigma[i], PP[i], BBphi[i], BBtheta[i], qq[i], qp[i]);
    else
      fprintf (file, "%e %e %e %e %e %e %e\n", rr[i], ssigma[i], -PP[i], BBphi[i], BBtheta[i], qq[i], qp[i]);
  
  fclose (file);
}

// ##################################################
// Function to calculate stepped pressure equilibrium
// ##################################################
void NonLinear::CalcStepped ()
{
  // Calculate pressures between control surfaces
  for (int i = 0; i <= Ncntr; i++)
    {
      PPc[i] = GetPQ (rrc[i]);
    }

  // Calculate parallel current profiles
  for (int i = 0; i <= Ngrid; i++)
    {
      sssigma[i] = GetSigma (rr[i]);
    }
  
  // Calculate stepped pressure profile
  int cntr = 1;
  for (int i = 0; i <= Ngrid; i++)
    {
      PPP[i] = PPc[cntr];

      if (rr[i] > rrc[cntr])
	cntr++;
    }
   
  // Calculate magnetic field and safety-factor profiles
  BBBphi[0]   = 1.;
  BBBtheta[0] = 0.;
  qqq[0]      = q0;

  double  r, h, t_err;
  int     rept; rhs_chooser = 2;
  double* y   = new double[2];
  double* err = new double[2];

  r     = eps;
  h     = h0;
  count = 0;

  y[0] = q0;
  y[1] = (epsa /q0) * r;

  cntr = 1;
  for (int i = 1; i <= Ngrid; i++)
    {
      while (rrc[cntr] < rr[i])
	{
	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rrc[cntr]);
	  RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
	  
	  // Apply jump condition at control surface
	  double q       = y[0];
	  double Btheta  = y[1];
	  double jump    = 2. * (PPc[cntr] - PPc[cntr+1]) * epsa*epsa * r*r /(epsa*epsa * r*r + q*q);
	  double Bt2p    = Btheta*Btheta + jump;
	  y[1]           = sqrt (Bt2p);
	  Btheta2j[cntr] = jump;
	  Bthetam[cntr]  = Btheta;
	  Bthetap[cntr]  = y[1];
	  cntr++;
	}
      
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      RK4RK5Fixed (2, r, y, err, rr[i] - r);
      
      qqq[i]      = y[0];
      BBBtheta[i] = y[1];
      BBBphi[i]   = qqq[i] * BBBtheta[i] /epsa /r;
    }

  delete[] y; delete[] err;

  // Output control surface data
  FILE* file = OpenFilew ("Output/Control.out");

  for (int i = 1; i < Ncntr; i++)
    fprintf (file, "%d %e %e %e %e %e\n", i, rrc[i], PPc[i], PPc[i+1], Bthetam[i], Bthetap[i]);
  
  fclose (file);
  
  // Output stepped pressure profiles
  file = OpenFilew ("Output/Stepped.out");

  for (int i = 0; i <= Ngrid; i++)
    fprintf (file, "%e %e %e %e %e %e\n", rr[i], sssigma[i], PPP[i], BBBphi[i], BBBtheta[i], qqq[i]);
  
  fclose (file);
}

// ###############################################################
// Function to calculate tearing eigenfunction (no pressure jumps)
// ###############################################################
void NonLinear::GetPsi (double _mpol, double _ntor, double rs, double& Delta, double* Psi, double* Psip, int verbose)
{
  // Set mode numbers
  mpol = _mpol;
  ntor = _ntor;
  
  // Launch solution from magnetic axis and integrate to rational surface
  double r, h, t_err, alpha, beta, gamma;
  int    rept; rhs_chooser = 3;
  double y[2], err[2];
  
  if (mpol > 1.e-6)
    {
      if (int(alphas) == 2)
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0 + mpol*mpol * nus /(mpol - ntor * q0)/(mpol - ntor * q0)) /mpol/mpol/mpol;
      else
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0                                                         ) /mpol/mpol/mpol;

      alpha = ((mpol + 2.) * ntor*ntor * epsa*epsa /mpol -             gamma) /((mpol+2.)*(mpol+2.) - 1.);
      beta  = (              ntor*ntor * epsa*epsa /mpol - (mpol+2.) * gamma) /((mpol+2.)*(mpol+2.) - 1.);
      
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r);
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r);
      h     = h0;
      count = 0;
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.);
      y[1]  = 2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.;
      h     = h0;
      count = 0;
    }

  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e                        psi = %11.4e  chi = %11.4e\n",
	    mpol, ntor, r, y[0], y[1]);

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);

  // Calculate coefficients of large and small solutions
  double H      = GetH (rs);
  double lambda = Getlambda (rs);
  double x      = - delta /rs;
  double a11    = 1. + lambda * x * (log (-x) - 1.);
  double a12    = x;
  double a21    = lambda * log (-x) /H;
  double a22    = 1. /H;
  double del    = a11 * a22 - a12 * a21;
  double Clm    = (  a22 * y[0] - a12 * y[1]) /del;
  double Csm    = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clm = %11.4e  Csm = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clm, Csm);

  // Launch solution from plasma boundary and integrate to rational surface
  double Ima  = gsl_sf_bessel_In (int (mpol),     epsa);
  double Kma  = gsl_sf_bessel_Kn (int (mpol),     epsa);
  double Im1a = gsl_sf_bessel_In (int (mpol) + 1, epsa);
  double Km1a = gsl_sf_bessel_Kn (int (mpol) + 1, epsa);
  double Impa =    Im1a + Ima * mpol /epsa;
  double Kmpa =  - Km1a + Kma * mpol /epsa;
  double Imb  = gsl_sf_bessel_In (int (mpol),     epsb);
  double Kmb  = gsl_sf_bessel_Kn (int (mpol),     epsb);
  double Im1b = gsl_sf_bessel_In (int (mpol) + 1, epsb);
  double Km1b = gsl_sf_bessel_Kn (int (mpol) + 1, epsb);
  double Impb =    Im1b + Imb * mpol /epsb;
  double Kmpb =  - Km1b + Kmb * mpol /epsb;
  
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb);
  y[1]  = Ima * Kmpb - Kma * Impb;
  h     = - h0;
  count = 0;

  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e                        psi = %11.4e  chi = %11.4e\n",
	    mpol, ntor, r, y[0], y[1]);

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + delta);
  RK4RK5Fixed (2, r, y, err, rs + delta - r);

  x          = delta /rs;
  a11        = 1. + lambda * x * (log (x) - 1.);
  a12        = x;
  a21        = lambda * log (x) /H;
  a22        = 1. /H;
  del        = a11 * a22 - a12 * a21;
  double Clp = (  a22 * y[0] - a12 * y[1]) /del;
  double Csp = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clp = %11.4e  Csp = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clp, Csp);

  // Calculate tearing stability index
  Delta = Csp /Clp - Csm /Clm;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  rs = %11.4e  Delta  = %11.4e\n",
	    mpol, ntor, rs, Delta);

  // -----------------------
  // Calculate eigenfunction
  // -----------------------
  double Norm = Clm, Norm0 = Clm, xi0 = - q0 /(q0 - 1.) /epsa;
  if (mpol > 1.e-6)
    {
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r) /Norm;
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r) /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  Psi[0]  = Norm0 * y[0] /xi0;
	  Psip[0] = Norm0 * y[1] /Getf(r) /xi0;
	}
      else
	{
	  Psi[0]  = y[0];
	  Psip[0] = y[1] /Getf(r);
	}
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.) /Norm;
      y[1]  = (2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.)         /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  Psi[0]  = Norm0 * y[0] /xi0;
	  Psip[0] = Norm0 * y[1] /Getf(r) /xi0;
	}
      else
	{
	  Psi[0]  = y[0];
	  Psip[0] = y[1] /Getf(r);
	}
    }

  for (int i = 1; i <= Ngrid; i++)
    {
      if (rr[i] < rs)
	{
	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      Psi[i]  = Norm0 * y[0] /xi0;
	      Psip[i] = Norm0 * y[1] /Getf(r) /xi0;
	    }
	  else
	    {
	      Psi[i]  = y[0];
	      Psip[i] = y[1] /Getf(r);
	    }
	}
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);
  
  x  = - delta /rs;
  a11 = 1. + lambda * x * (log (-x) - 1.);
  a12 = x;
  a21 = lambda * log (-x) /H;
  a22 = 1. /H;
  del = a11 * a22 - a12 * a21;
  Clm = (  a22 * y[0] - a12 * y[1]) /del;
  Csm = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clm = %11.4e  Csm = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clm, Csm);

  Norm  = Clp;
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb) /Norm;
  y[1]  = (Ima * Kmpb - Kma * Impb)                 /Norm;
  h     = - h0;
  count = 0;

  if (int(mpol) == 1 && int(ntor) == 1)
    {
      Psi[Ngrid]  = Norm0 * y[0] /xi0;
      Psip[Ngrid] = Norm0 * y[1] /Getf(r) /xi0;
    }
  else
    {
      Psi[Ngrid]  = y[0];
      Psip[Ngrid] = y[1] /Getf(r);
    }

  int cntr = 1;
  for (int i = Ngrid-1; i >= 0; i--)
    {
      if (rr[i] > rs)
	{
	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r > rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      Psi[i]  = Norm0 * y[0] /xi0;
	      Psip[i] = Norm0 * y[1] /Getf(r) /xi0;
	    }
	  else
	    {
	      Psi[i]  = y[0];
	      Psip[i] = y[1] /Getf(r);
	    }
	}
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + delta);
  RK4RK5Fixed (2, r, y, err, rs + delta - r);

  x   = delta /rs;
  a11 = 1. + lambda * x * (log (x) - 1.);
  a12 = x;
  a21 = lambda * log (x) /H;
  a22 = 1. /H;
  del = a11 * a22 - a12 * a21;
  Clp = (  a22 * y[0] - a12 * y[1]) /del;
  Csp = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clp = %11.4e  Csp = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clp, Csp);
}

// ###################################################################
// Function to calculate tearing stability index (with pressure jumps)
// ###################################################################
void NonLinear::GetDelta (double _mpol, double _ntor, double rs, double* JJ, double& Delta, double& Psirv)
{
  // Set mode numbers
  mpol = _mpol;
  ntor = _ntor;
  
  // Launch solution from magnetic axis and integrate to rational surface
  double r, h, t_err, alpha, beta, gamma;
  int    rept; rhs_chooser = 3;
  double y[2], err[2];

  if (mpol > 1.e-6)
    {
      if (int(alphas) == 2)
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0 + mpol*mpol * nus /(mpol - ntor * q0)/(mpol - ntor * q0)) /mpol/mpol/mpol;
      else
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0                                                         ) /mpol/mpol/mpol;

      alpha = ((mpol + 2.) * ntor*ntor * epsa*epsa /mpol -             gamma) /((mpol+2.)*(mpol+2.) - 1.);
      beta  = (              ntor*ntor * epsa*epsa /mpol - (mpol+2.) * gamma) /((mpol+2.)*(mpol+2.) - 1.);
      
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r);
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r);
      h     = h0;
      count = 0;
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.);
      y[1]  = 2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.;
      h     = h0;
      count = 0;
    }

  int cntr = 1;
  while (rrc[cntr] < rs - delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      // Apply jump conditions at control surface
      y[1] += JJ[cntr] * y[0];
      cntr++;
    }
  
  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);

  // Calculate coefficients of large and small solutions
  double H      = GetH (rs);
  double lambda = Getlambda (rs);
  double x      = - delta /rs;
  double a11    = 1. + lambda * x * (log (-x) - 1.);
  double a12    = x;
  double a21    = lambda * log (-x) /H;
  double a22    = 1. /H;
  double del    = a11 * a22 - a12 * a21;
  double Clm    = (  a22 * y[0] - a12 * y[1]) /del;
  double Csm    = (- a21 * y[0] + a11 * y[1]) /del;

  // Launch solution from plasma boundary and integrate to rational surface
  double Ima  = gsl_sf_bessel_In (int (mpol),     epsa);
  double Kma  = gsl_sf_bessel_Kn (int (mpol),     epsa);
  double Im1a = gsl_sf_bessel_In (int (mpol) + 1, epsa);
  double Km1a = gsl_sf_bessel_Kn (int (mpol) + 1, epsa);
  double Impa =    Im1a + Ima * mpol /epsa;
  double Kmpa =  - Km1a + Kma * mpol /epsa;
  double Imb  = gsl_sf_bessel_In (int (mpol),     epsb);
  double Kmb  = gsl_sf_bessel_Kn (int (mpol),     epsb);
  double Im1b = gsl_sf_bessel_In (int (mpol) + 1, epsb);
  double Km1b = gsl_sf_bessel_Kn (int (mpol) + 1, epsb);
  double Impb =    Im1b + Imb * mpol /epsb;
  double Kmpb =  - Km1b + Kmb * mpol /epsb;
  
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb);
  y[1]  = Ima * Kmpb - Kma * Impb;
  h     = - h0;
  count = 0;

  cntr = Ncntr - 1;
  while (rrc[cntr] > rs + delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      // Apply jump conditions at control surface
      y[1] -= JJ[cntr] * y[0];
      cntr--;
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + delta);
  RK4RK5Fixed (2, r, y, err, rs + delta - r);
  
  x          = delta /rs;
  a11        = 1. + lambda * x * (log (x) - 1.);
  a12        = x;
  a21        = lambda * log (x) /H;
  a22        = 1. /H;
  del        = a11 * a22 - a12 * a21;
  double Clp = (  a22 * y[0] - a12 * y[1]) /del;
  double Csp = (- a21 * y[0] + a11 * y[1]) /del;
 
  // Calculate tearing stability index
  Delta = Csp /Clp - Csm /Clm;

  // -----------------------
  // Calculate eigenfunction
  // -----------------------
  double* PPsi  = new double[Ngrid+1];
  double Norm = Clm, Norm0 = Clm, F, Fp = GetFp (rs), xi0 = - q0 /(q0 - 1.) /epsa;
  if (mpol > 1.e-6)
    {
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r) /Norm;
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r) /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  PPsi[0]  = Norm0 * y[0] /xi0;
	}
      else
	{
	  PPsi[0]  = y[0];
	}
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.) /Norm;
      y[1]  = (2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.)         /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  PPsi[0]  = Norm0 * y[0] /xi0;
	}
      else
	{
	  PPsi[0]  = y[0];
	}
    }

  cntr = 1;
  for (int i = 1; i <= Ngrid; i++)
    {
      if (rr[i] < rs)
	{
	  while (rrc[cntr] < rr[i])
	    {
	      do
		{
		  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
		}
	      while (r < rrc[cntr]);
	      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
	      y[1] += JJ[cntr] * y[0];
	      cntr++;
	    }

	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      PPsi[i]  = Norm0 * y[0] /xi0;
	    }
	  else
	    {
	      PPsi[i]  = y[0];
	    }
	}
    }

  while (rrc[cntr] < rs - delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      y[1] += JJ[cntr] * y[0];
      cntr++;
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);
  
  x   = - delta /rs;
  a11 = 1. + lambda * x * (log (-x) - 1.);
  a12 = x;
  a21 = lambda * log (-x) /H;
  a22 = 1. /H;
  del = a11 * a22 - a12 * a21;
  Clm = (  a22 * y[0] - a12 * y[1]) /del;
  Csm = (- a21 * y[0] + a11 * y[1]) /del;

  Norm  = Clp;
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb) /Norm;
  y[1]  = (Ima * Kmpb - Kma * Impb)                 /Norm;
  h     = - h0;
  count = 0;
  if (int(mpol) == 1 && int(ntor) == 1)
    {
      PPsi[Ngrid]  = Norm0 * y[0] /xi0;
    }
  else
    {
      PPsi[Ngrid]  = y[0];
    }

  cntr = Ncntr - 1;
  for (int i = Ngrid-1; i >= 0; i--)
    {
      if (rr[i] > rs)
	{
	  while (rrc[cntr] > rr[i])
	    {
	      do
		{
		  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
		}
	      while (r > rrc[cntr]);
	      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);

	      y[1] -= JJ[cntr] * y[0];
	      cntr--;
	    }
	  
	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r > rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      PPsi[i]  = Norm0 * y[0] /xi0;
	    }
	  else
	    {
	      PPsi[i]  = y[0];
	    }
	}
    }

  // Find reversal surface
  double rv = Findrs (0.);

  // Get eigenfunction value at reversal surface
  int i;
  for (i = 0; i <= Ngrid; i++)
    if (rr[i] <= rv && rr[i+1] > rv)
      break;
  Psirv = (PPsi[i] * (rr[i+1] - rv) + PPsi[i+1] * (rv - rr[i])) /(rr[i+1] - rr[i]);

  delete[] PPsi; 
 }

// #################################################################
// Function to calculate tearing eigenfunction (with pressure jumps)
// #################################################################
void NonLinear::GetPsiStepped (double _mpol, double _ntor, double rs, double W, double* JJ, double& Delta, double* PPsi, double* PPsip, double* Xi, int verbose)
{
  // Set mode numbers
  mpol = _mpol;
  ntor = _ntor;
  
  // Launch solution from magnetic axis and integrate to rational surface
  double r, h, t_err, alpha, beta, gamma;
  int    rept; rhs_chooser = 3;
  double y[2], err[2];

  if (mpol > 1.e-6)
    {
      if (int(alphas) == 2)
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0 + mpol*mpol * nus /(mpol - ntor * q0)/(mpol - ntor * q0)) /mpol/mpol/mpol;
      else
	gamma = 4. * (mpol - ntor * q0) * (epsa*epsa /q0/q0                                                         ) /mpol/mpol/mpol;

      alpha = ((mpol + 2.) * ntor*ntor * epsa*epsa /mpol -             gamma) /((mpol+2.)*(mpol+2.) - 1.);
      beta  = (              ntor*ntor * epsa*epsa /mpol - (mpol+2.) * gamma) /((mpol+2.)*(mpol+2.) - 1.);
      
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r);
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r);
      h     = h0;
      count = 0;
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.);
      y[1]  = 2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.;
      h     = h0;
      count = 0;
    }

  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e                        psi = %11.4e  chi = %11.4e\n",
	    mpol, ntor, r, y[0], y[1]);

  int cntr = 1;
  while (rrc[cntr] < rs - delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      // Apply jump conditions at control surface
      y[1] += JJ[cntr] * y[0];
      cntr++;
    }
  
  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);

  // Calculate coefficients of large and small solutions
  double H      = GetH (rs);
  double lambda = Getlambda (rs);
  double x      = - delta /rs;
  double a11    = 1. + lambda * x * (log (-x) - 1.);
  double a12    = x;
  double a21    = lambda * log (-x) /H;
  double a22    = 1. /H;
  double del    = a11 * a22 - a12 * a21;
  double Clm    = (  a22 * y[0] - a12 * y[1]) /del;
  double Csm    = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clm = %11.4e  Csm = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clm, Csm);

  // Launch solution from plasma boundary and integrate to rational surface
  double Ima  = gsl_sf_bessel_In (int (mpol),     epsa);
  double Kma  = gsl_sf_bessel_Kn (int (mpol),     epsa);
  double Im1a = gsl_sf_bessel_In (int (mpol) + 1, epsa);
  double Km1a = gsl_sf_bessel_Kn (int (mpol) + 1, epsa);
  double Impa =    Im1a + Ima * mpol /epsa;
  double Kmpa =  - Km1a + Kma * mpol /epsa;
  double Imb  = gsl_sf_bessel_In (int (mpol),     epsb);
  double Kmb  = gsl_sf_bessel_Kn (int (mpol),     epsb);
  double Im1b = gsl_sf_bessel_In (int (mpol) + 1, epsb);
  double Km1b = gsl_sf_bessel_Kn (int (mpol) + 1, epsb);
  double Impb =    Im1b + Imb * mpol /epsb;
  double Kmpb =  - Km1b + Kmb * mpol /epsb;
  
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb);
  y[1]  = Ima * Kmpb - Kma * Impb;
  h     = - h0;
  count = 0;

  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e                        psi = %11.4e  chi = %11.4e\n",
	    mpol, ntor, r, y[0], y[1]);

  cntr = Ncntr - 1;
  while (rrc[cntr] > rs + delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      // Apply jump conditions at control surface
      y[1] -= JJ[cntr] * y[0];
      cntr--;
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + delta);
  RK4RK5Fixed (2, r, y, err, rs + delta - r);
  
  x          = delta /rs;
  a11        = 1. + lambda * x * (log (x) - 1.);
  a12        = x;
  a21        = lambda * log (x) /H;
  a22        = 1. /H;
  del        = a11 * a22 - a12 * a21;
  double Clp = (  a22 * y[0] - a12 * y[1]) /del;
  double Csp = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clp = %11.4e  Csp = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clp, Csp);

  // Calculate tearing stability index
  Delta = Csp /Clp - Csm /Clm;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  rs = %11.4e  Delta  = %11.4e\n",
	    mpol, ntor, rs, Delta);

  // -----------------------
  // Calculate eigenfunction
  // -----------------------
  double Norm = Clm, Norm0 = Clm, F, Fp = GetFp (rs), xi0 = - q0 /(q0 - 1.) /epsa;
  if (mpol > 1.e-6)
    {
      r     = eps;
      y[0]  = pow (r, mpol) * (1.       + alpha *r*r) /Norm;
      y[1]  = pow (r, mpol) * (1. /mpol + beta * r*r) /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  PPsi[0]  = Norm0 * y[0] /xi0;
	  PPsip[0] = Norm0 * y[1] /Getf(r) /xi0;
	}
      else
	{
	  PPsi[0]  = y[0];
	  PPsip[0] = y[1] /Getf(r);
	}
      F     = Getq (r) /r /epsa /GetBphi (r);
      Xi[0] = PPsi[0] * F * Getimnq (r, rs, W);
    }
  else
    {
      r     = eps;
      y[0]  = r*r * (1. + ntor*ntor * epsa*epsa * (1. - 4. /ntor/ntor /q0/q0) * r*r /8.) /Norm;
      y[1]  = (2. /ntor/ntor /epsa/epsa + (1. - 4. /ntor/ntor /q0/q0) * r*r /2.)         /Norm;
      h     = h0;
      count = 0;
      if (int(mpol) == 1 && int(ntor) == 1)
	{
	  PPsi[0]  = Norm0 * y[0] /xi0;
	  PPsip[0] = Norm0 * y[1] /Getf(r) /xi0;
	}
      else
	{
	  PPsi[0]  = y[0];
	  PPsip[0] = y[1] /Getf(r);
	}
      F     = Getq (r) /r /epsa /GetBphi (r);
      Xi[0] = PPsi[0] * F * Getimnq (r, rs, W);
    }

  cntr = 1;
  for (int i = 1; i <= Ngrid; i++)
    {
      if (rr[i] < rs)
	{
	  while (rrc[cntr] < rr[i])
	    {
	      do
		{
		  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
		}
	      while (r < rrc[cntr]);
	      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
	      y[1] += JJ[cntr] * y[0];
	      cntr++;
	    }

	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      PPsi[i]  = Norm0 * y[0] /xi0;
	      PPsip[i] = Norm0 * y[1] /Getf(r) /xi0;
	    }
	  else
	    {
	      PPsi[i]  = y[0];
	      PPsip[i] = y[1] /Getf(r);
	    }
	  F     = Getq (r) /r /epsa /GetBphi (r);
	  Xi[i] = PPsi[i] * F * Getimnq (r, rs, W);
	}
    }

  while (rrc[cntr] < rs - delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      y[1] += JJ[cntr] * y[0];
      cntr++;
    }

  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - delta);
  RK4RK5Fixed (2, r, y, err, rs - delta - r);
  
  x   = - delta /rs;
  a11 = 1. + lambda * x * (log (-x) - 1.);
  a12 = x;
  a21 = lambda * log (-x) /H;
  a22 = 1. /H;
  del = a11 * a22 - a12 * a21;
  Clm = (  a22 * y[0] - a12 * y[1]) /del;
  Csm = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clm = %11.4e  Csm = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clm, Csm);

  Norm  = Clp;
  r     = 1.;
  y[0]  = ntor * epsa * (Impa * Kmpb - Kmpa * Impb) /Norm;
  y[1]  = (Ima * Kmpb - Kma * Impb)                 /Norm;
  h     = - h0;
  count = 0;
  if (int(mpol) == 1 && int(ntor) == 1)
    {
      PPsi[Ngrid]  = Norm0 * y[0] /xi0;
      PPsip[Ngrid] = Norm0 * y[1] /Getf(r) /xi0;
    }
  else
    {
      PPsi[Ngrid]  = y[0];
      PPsip[Ngrid] = y[1] /Getf(r);
    }
  F         = Getq (r) /r /epsa /GetBphi (r);
  Xi[Ngrid] = PPsi[Ngrid] * F * Getimnq (r, rs, W);

  cntr = Ncntr - 1;
  for (int i = Ngrid-1; i >= 0; i--)
    {
      if (rr[i] > rs)
	{
	  while (rrc[cntr] > rr[i])
	    {
	      do
		{
		  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
		}
	      while (r > rrc[cntr]);
	      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);

	      y[1] -= JJ[cntr] * y[0];
	      cntr--;
	    }
	  
	  do
	    {
	      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r > rr[i]);
	  RK4RK5Fixed (2, r, y, err, rr[i] - r);

	  if (int(mpol) == 1 && int(ntor) == 1)
	    {
	      PPsi[i]  = Norm0 * y[0] /xi0;
	      PPsip[i] = Norm0 * y[1] /Getf(r) /xi0;
	    }
	  else
	    {
	      PPsi[i]  = y[0];
	      PPsip[i] = y[1] /Getf(r);
	    }
	  F     = Getq (r) /r /epsa /GetBphi (r);
	  Xi[i] = PPsi[i] * F * Getimnq (r, rs, W);
	}
    }

   while (rrc[cntr] > rs + delta)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rrc[cntr]);
      RK4RK5Fixed (2, r, y, err, rrc[cntr] - r);
      
      y[1] -= JJ[cntr] * y[0];
      cntr--;
    }
  
  do
    {
      RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + delta);
  RK4RK5Fixed (2, r, y, err, rs + delta - r);

  x   = delta /rs;
  a11 = 1. + lambda * x * (log (x) - 1.);
  a12 = x;
  a21 = lambda * log (x) /H;
  a22 = 1. /H;
  del = a11 * a22 - a12 * a21;
  Clp = (  a22 * y[0] - a12 * y[1]) /del;
  Csp = (- a21 * y[0] + a11 * y[1]) /del;
  if (verbose)
    printf ("             mpol = %11.4e  ntor = %11.4e  r  = %11.4e  r - rs = %11.4e  psi = %11.4e  chi = %11.4e  Clp = %11.4e  Csp = %11.4e\n",
	    mpol, ntor, rs, r - rs, y[0], y[1], Clp, Csp);
}

// ######################################
// Function to calculate overlap function
// ######################################
void NonLinear::GetTau ()
{
  Tau[0]  = 0.;
  Taup[0] = 0.;
  
  for (int i = 1; i <= Ngrid; i++)
    {
      double r = rr[i];

      double F1  = GetF    (r, mpol1, ntor1);
      double F2  = GetF    (r, mpol2, ntor2);
      double F3  = GetF    (r, mpol3, ntor3);
      double iF1 = GetinvF (r, rs1, W1, mpol1, ntor1);
      double iF2 = GetinvF (r, rs2, W2, mpol2, ntor2);
      double iF3 = GetinvF (r, rs3, W3, mpol3, ntor3);
      double G1  = GetG    (r, mpol1, ntor1);
      double G2  = GetG    (r, mpol2, ntor2);
      double G3  = GetG    (r, mpol3, ntor3);
      double H1  = GetH    (r, mpol1, ntor1);
      double H2  = GetH    (r, mpol2, ntor2);
      double H3  = GetH    (r, mpol3, ntor3);

      double sigma  = GetSigma  (r);
      double sigmap = GetSigmap (r);
      double Btheta = GetBtheta (r);
      double Bphi   = GetBphi   (r);

      double f1;
      f1  = r * PPsip1[i] * PPsi2[i] * PPsi3[i] * G1 * iF2 * iF3 /H1;
      f1 += r * PPsip2[i] * PPsi3[i] * PPsi1[i] * G2 * iF3 * iF1 /H2;
      f1 += r * PPsip3[i] * PPsi1[i] * PPsi2[i] * G3 * iF1 * iF2 /H3;

      double f2;
      f2  = F1 * iF2 * iF3 /H1;
      f2 += F2 * iF3 * iF1 /H2;
      f2 += F3 * iF1 * iF2 /H3;

      f2 *= r * sigma * PPsi1[i] * PPsi2[i] * PPsi3[i]; 

      double f3;
      f3 = (2. * Btheta * Bphi - r * sigma * (Btheta*Btheta + Bphi*Bphi)) * iF1 * iF2 * iF3;

      f3 *= PPsi1[i] * PPsi2[i] * PPsi3[i]; 

      double f4 = sigmap * (f1 + f2 + f3);

      Taup[i] = f4;

      double f10;
      f10  = r * Psip1[i] * Psi2[i] * Psi3[i] * G1 * iF2 * iF3 /H1;
      f10 += r * Psip2[i] * Psi3[i] * Psi1[i] * G2 * iF3 * iF1 /H2;
      f10 += r * Psip3[i] * Psi1[i] * Psi2[i] * G3 * iF1 * iF2 /H3;

      double f20;
      f20  = F1 * iF2 * iF3 /H1;
      f20 += F2 * iF3 * iF1 /H2;
      f20 += F3 * iF1 * iF2 /H3;

      f20 *= r * sigma * Psi1[i] * Psi2[i] * Psi3[i]; 

      double f30;
      f30 = (2. * Btheta * Bphi - r * sigma * (Btheta*Btheta + Bphi*Bphi)) * iF1 * iF2 * iF3;

      f30 *= Psi1[i] * Psi2[i] * Psi3[i]; 

      double f40 = sigmap * (f10 + f20 + f30);

      Tau[i] = f40;
    }
}

// ######################################
// Function to calculate overlap integral
// ######################################
void NonLinear::CalcOverlap ()
{
  double r, h, t_err;
  int    rept; rhs_chooser = 4;
  double y[2], err[2];

  r     = eps;
  h     = h0;
  count = 0;

  y[0] = 0.;
  y[1] = 0.;

  for (int i = 1; i <= Ngrid; i++)
    {
      do
	{
	  RK4RK5Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      RK4RK5Fixed (2, r, y, err, rr[i] - r);

      Over[i]  = y[0];
      Overp[i] = y[1];
    }
}

// ##########q#############################
// Function to find rational surface radius
// ########################################
double NonLinear::Findrs (double qs)
{
  double rs = -1.;
  
  int N = 1000;
  for (int i = 0; i < N; i++)
    {
      double r1 = double (i)   /double (N);
      double r2 = double (i+1) /double (N);
      double q1 = Getq (r1);
      double q2 = Getq (r2);
      //printf ("i = %3d  r1 = %11.4e  q1 = %11.4e  r2 = %11.4e  q2 = %11.4e\n", i, r1, q1, r2, q2);

      if ((q1 - qs) * (q2 - qs) < 0.)
	{
	  double rn = r1, qn;
	  int count = 0;
	  do
	    {
	      qn = Getq (rn);
	      double qpn = Getqp (rn);
	      
	      //printf ("rn = %11.4e  qn = %11.4e  residual = %11.4e\n", rn, qn, fabs (qn - qs));
	      
	      rn += (qs - qn) /qpn;
	      count++;
	    }
	  while (fabs (qn - qs) > 1.e-15 && count < 100);
	  
	  if (isnan (rn))
	    {
	      printf ("NonLinear:Findrs - Error: rn is NaN\n");
	      exit (1);
	    }
	  else
	    rs = rn;
	}
    }

  return rs;
}

// #############################
// Function to return value of F
// #############################
double NonLinear::GetF (double r, double _mpol, double _ntor)
{
  double Bphi   = GetBphi (r);
  double Btheta = GetBtheta (r);

  return _mpol * Btheta - _ntor * epsa * r * Bphi;
}

// ###########################################
// Function to return regularized value of 1/F
// ###########################################
double NonLinear::GetinvF (double r, double rs, double W, double _mpol, double _ntor)
{
  double q    = Getq (r);
  double qp   = Getqp (rs);
  double Bphi = GetBphi (r);
  double mnq  = _mpol - _ntor * q;

  return (q /epsa /r /Bphi) * mnq /(mnq*mnq + _ntor*_ntor * qp*qp * W*W);
}

// #############################
// Function to return value of G
// #############################
double NonLinear::GetG (double r, double _mpol, double _ntor)
{
  double Bphi   = GetBphi (r);
  double Btheta = GetBtheta (r);

  return _ntor * epsa * r * Btheta + _mpol * Bphi;
}

// #############################
// Function to return value of H
// #############################
double NonLinear::GetH (double r, double _mpol, double _ntor)
{
  return _mpol * _mpol + _ntor * _ntor * epsa * epsa * r * r;
}

// ##########################################################
// Function to return the values of J at the control surfaces
// ##########################################################
void NonLinear::GetJ (double _mpol, double _ntor, double rs, double W, double* JJ)
{
  mpol = _mpol;
  ntor = _ntor;
  
  for (int i = 0; i < Ncntr; i++)
    {
      double K = GetK (rrc[i], rs, W);
      double H = GetH (rrc[i]);
      JJ[i]    = K * Btheta2j[i] /Bthetam[i] /Bthetap[i] /H;
    }
}

// #############################
// Function to return value of K
// #############################
double NonLinear::GetK (double r, double rs, double W)
{
  double q     = Getq (r);
  double sigma = GetSigma (r);
  double H     = GetH (r);
  double imnq  = Getimnq (r, rs, W);
  double mnh   = mpol * ntor * q + ntor*ntor * epsa*epsa * r*r;

  return imnq * (mpol + mnh * (sigma /2. /ntor /epsa - mpol /H)) - imnq*imnq * mnh;
}

// ########################################
// Function to return value of 1 /(m - n q)
// ########################################
double NonLinear::Getimnq (double r, double rs, double W)
{
  double q   = Getq (r);
  double qp  = Getqp (rs);
  double mnq = mpol - ntor * q;

  return mnq /(mnq * mnq + ntor*ntor * qp*qp * W*W);
}

// #############################
// Function to return value of f
// #############################
double NonLinear::Getf (double r)
{
  double H = GetH (r);

  return r /H;
}

// #############################
// Function to return value of g
// #############################
double NonLinear::Getg (double r)
{
  double F      = GetF (r);
  double G      = GetG (r);
  double H      = GetH (r);
  double sigma  = GetSigma (r);
  double sigmap = GetSigmap (r);

  return 1. /r + 2. * mpol * ntor * epsa * r * sigma /H/H - r * sigma*sigma /H + G * r * sigmap /H /F;
}

// ##################################
// Function to return value of lambda
// ##################################
double NonLinear::Getlambda (double r)
{
  double G      = GetG (r);
  double Fp     = GetFp (r);
  double sigmap = GetSigmap (r);

  return G * r * sigmap /Fp;
}

// #############################
// Function to return value of F
// #############################
double NonLinear::GetF (double r)
{
  double Bphi   = GetBphi (r);
  double Btheta = GetBtheta (r);

  return mpol * Btheta - ntor * epsa * r * Bphi;
}

// ##############################
// Function to return value of Fp
// ##############################
double NonLinear::GetFp (double r)
{
  double Bphi   = GetBphi (r);
  double Btheta = GetBtheta (r);
  double sigma  = GetSigma (r);
  double G      = GetG (r);

  return sigma * G - mpol * Btheta /r - ntor * epsa * Bphi;
}

// #############################
// Function to return value of G
// #############################
double NonLinear::GetG (double r)
{
  double Bphi   = GetBphi (r);
  double Btheta = GetBtheta (r);

  return ntor * epsa * r * Btheta + mpol * Bphi;
}

// #############################
// Function to return value of H
// #############################
double NonLinear::GetH (double r)
{
  return mpol * mpol + ntor * ntor * epsa * epsa * r * r;
}

// ###################################################
// Function to return value of toroidal magnetic field
// ###################################################
double NonLinear::GetBphi (double r)
{
  double dr = 1. /double (Ngrid);
  int    n  = int (r /dr);

  if (n <= 0)
    n = 1;
  if (n >= Ngrid)
    n = Ngrid - 1;

  double s = (r - rr[n]) /dr; 
  
  return 0.5*s*(s-1.) * BBphi[n-1] - (s+1.)*(s-1.) * BBphi[n] + 0.5*s*(s+1.) * BBphi[n+1];
}

// ###################################################
// Function to return value of poloidal magnetic field
// ###################################################
double NonLinear::GetBtheta (double r)
{
  double dr = 1. /double (Ngrid);
  int    n  = int (r /dr);

  if (n <= 0)
    n = 1;
  if (n >= Ngrid)
    n = Ngrid - 1;

  double s = (r - rr[n]) /dr; 
  
  return 0.5*s*(s-1.) * BBtheta[n-1] - (s+1.)*(s-1.) * BBtheta[n] + 0.5*s*(s+1.) * BBtheta[n+1];
}

// #########################################
// Function to return value of safety-factor
// #########################################
double NonLinear::Getq (double r)
{
  double dr = 1. /double (Ngrid);
  int    n  = int (r /dr);

  if (n <= 0)
    n = 1;
  if (n >= Ngrid)
    n = Ngrid - 1;

  double s = (r - rr[n]) /dr; 
  
  return 0.5*s*(s-1.) * qq[n-1] - (s+1.)*(s-1.) * qq[n] + 0.5*s*(s+1.) * qq[n+1];
}

// ##################################################
// Function to return value of safety-factor gradient
// ##################################################
double NonLinear::Getqp (double r)
{
  double dr = 1. /double (Ngrid);
  int    n  = int (r /dr);

  if (n <= 0)
    n = 1;
  if (n >= Ngrid)
    n = Ngrid - 1;

  double s = (r - rr[n]) /dr; 

  return 0.5*s*(s-1.) * qp[n-1] - (s+1.)*(s-1.) * qp[n] + 0.5*s*(s+1.) * qp[n+1];
}

// ############################################
// Function to get equilibrium parallel current
// ############################################
double NonLinear::GetSigma (double r)
{
  if (r >= 1.)
    return 0;
  else
    return (2. * epsa /q0) * pow (1. - pow (r, alphas), nus);
}

// #####################################################
// Function to get equilibrium parallel current gradient
// #####################################################
double NonLinear::GetSigmap (double r)
{
  if (r >= 1.)
    return 0;
  else
    return - (2. * epsa /q0) * alphas * nus * pow (r, alphas - 1.) * pow (1. - pow (r, alphas), nus - 1.);
}

// ####################################
// Function to get equilibrium pressure
// ####################################
double NonLinear::GetP (double r)
{
  if (r >= 1.)
    return 0;
  else
    return (beta0 /2.) * pow (1. - pow (r, alphap), nup);
}

// #############################################
// Function to get equilibrium pressure gradient
// #############################################
double NonLinear::GetPp (double r)
{
  if (r >= 1.)
    return 0;
  else
    return - (beta0 /2.) * alphap * nup * pow (r, alphap - 1.) * pow (1. - pow (r, alphap), nup - 1.);
}

// #############################################
// Function to get modified equilibrium pressure
// #############################################
double NonLinear::GetPQ (double r)
{
  if (TorCur && r < 1.)
    {
      return gsl_spline_eval (SPQ, r, APQ);
    }
  else
    return GetP (r);
}

// ######################################################
// Function to get modified equilibrium pressure gradient
// ######################################################
double NonLinear::GetPQp (double r)
{
  if (TorCur)
    {
      double q = Getq (r);
      return (1. - q*q) * GetPp (r);
    }
  else
    return GetPp (r);
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void NonLinear::Rhs (double r, double* y, double* dydr)
{
  if (rhs_chooser == 0)
    {
      double sigma  = GetSigma (r);
      double Pp     = GetPp (r);
      double Bphi   = y[0];
      double Btheta = y[1];
      double q      = epsa * r * Bphi /Btheta;
      
      dydr[0] =             - sigma * Btheta;
      dydr[1] = - Btheta /r + sigma * Bphi;
      dydr[2] = (1. - q*q) * Pp;
      dydr[3] = Bphi * r;
    }
  else if (rhs_chooser == 1)
    {
      double sigma  = GetSigma (r);
      double Pp     = GetPQp (r);
      double Bphi   = y[0];
      double Btheta = y[1];
      double q      = epsa * r * Bphi /Btheta;

      dydr[0] =             - sigma * Btheta - Pp * Bphi   /(Btheta*Btheta + Bphi*Bphi);
      dydr[1] = - Btheta /r + sigma * Bphi   - Pp * Btheta /(Btheta*Btheta + Bphi*Bphi);
      dydr[2] = Bphi * r;
    }
  else if (rhs_chooser == 2)
    {
      double sigma  = GetSigma (r);
      double q      = y[0];
      double Btheta = y[1];
      
      dydr[0] = (2. - sigma * q /epsa) * q /r - epsa * r * sigma;
      dydr[1] = (sigma * q /epsa - 1.) * Btheta /r;
    }
  else if (rhs_chooser == 3)
    {
      double f = Getf (r);
      double g = Getg (r);

      dydr[0] = y[1] /f;
      dydr[1] = g * y[0];
    }
  else if (rhs_chooser == 4)
    {
      dydr[0] = gsl_spline_eval (Stau,  r, Atau);
      dydr[1] = gsl_spline_eval (Staup, r, Ataup);
    }
}

// #######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step length fourth-order/fifth-order Runge-Kutta scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... safety factor
//     T       ... step length cannot change by more than this factor from step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step length
//     h_max   ... maximum allowable step length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step while attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
// #######################################################################
void NonLinear::RK4RK5Adaptive (int neqns, double& x, double* y, double& h, 
			     double& t_err, double acc, double S, double T, int& rept,
			     int maxrept, double h_min, double h_max, int flag, 
			     int diag, FILE* file)
{
  double* y0  = new double[neqns];
  double* Err = new double[neqns];
  double  hin = h;

  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0[i] = y[i];

  // Take RK4/RK5 step 
  RK4RK5Fixed (neqns, x, y, Err, h);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err[i] / y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs (Err[i] / y[i]);
	  err2  = fabs (Err[i]);
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc / t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc / t_err), 0.25);

  // Prevent step length from changing by more than factor T
  if (h_est / h > T)
    h *= T;
  else if (h_est / h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Prevent step length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %11.4e hin = %11.4e err = %11.4e acc = %11.4e hout = %11.4e count = %3d\n", 
	     x, hin, t_err, acc, h, count);


  // Check if truncation error acceptable
  if ((t_err <= acc) || (count >= maxrept))
    {
      // If truncation error acceptable take step 
      rept  = count;
      count = 0;
    }
  else 
    {
      // If truncation error unacceptable repeat step 
      count++;
      x = x0;
      for (int i = 0; i < neqns; i++)
	y[i] = y0[i];
      RK4RK5Adaptive (neqns, x, y, h, t_err, acc, S, T, rept, 
		      maxrept, h_min, h_max, flag, diag, file);
    }

  delete[] y0; 
  delete[] Err;
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step length fourth-order/fifth-order Runge-Kutta scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables 
//     err     ... array of errors
//     h       ... step length
//     
// #####################################################################
void NonLinear::RK4RK5Fixed (int neqns, double& x, double* y, double* err, double h)
{
  double* dydx = new double[neqns];
  double* k1   = new double[neqns];
  double* k2   = new double[neqns];
  double* k3   = new double[neqns];
  double* k4   = new double[neqns];
  double* k5   = new double[neqns];
  double* k6   = new double[neqns];
  double* f    = new double[neqns];

  // First stage
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  Rhs (x + aa6 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k6[i] = h * dydx[i];
    }

  // Actual step 
  for (int i = 0; i < neqns; i++)
    {
      y  [i] = y[i] + cc1 * k1[i] + cc3 * k3[i] + cc4 * k4[i]               + cc6 * k6[i];
      err[i] =        ca1 * k1[i] + ca3 * k3[i] + ca4 * k4[i] + ca5 * k5[i] + ca6 * k6[i];
    }
  x += h;

  delete[] dydx; delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6; delete[] f;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* NonLinear::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

