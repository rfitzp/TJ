// Equilibrium.cpp

#include "LightEquilibrium.h"
#include "Equilibrium.h"

// ###########
// Constructor
// ###########
Equilibrium::Equilibrium ()
{
  // --------------------------------------------------
  // Ensure that directory ../Outputs/Equilibrium exits
  // --------------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Equilibrium"))
    {
      exit (1);
    }
  
  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/Equilibrium.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  SRC  = JSONData["SRC"] .get<int>    ();
  qc   = JSONData["qc"]  .get<double> ();
  qa   = JSONData["qa"]  .get<double> ();
  pc   = JSONData["pc"]  .get<double> ();
  mu   = JSONData["mu"]  .get<double> ();
  epsa = JSONData["epsa"].get<double> ();
  eps  = JSONData["eps"] .get<double> ();
  Ns   = JSONData["Ns"]  .get<int>    ();
  Nr   = JSONData["Nr"]  .get<int>    ();
  Nf   = JSONData["Nf"]  .get<int>    ();
  Nw   = JSONData["Nw"]  .get<int>    ();
  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();
  tilt = JSONData["tilt"].get<double> ();

  for (const auto& number : JSONData["Hna"])
    {
      Hna.push_back (number.get<double> ());
    }
  for (const auto& number : JSONData["Vna"])
    {
      Vna.push_back (number.get<double> ());
    }

  B0    = JSONData["B0"]   .get<double> ();
  R0    = JSONData["R0"]   .get<double> ();
  n0    = JSONData["n0"]   .get<double> ();
  alpha = JSONData["alpha"].get<double> ();
  Teped = JSONData["Teped"].get<double> ();
  neped = JSONData["neped"].get<double> ();
  
  JSONFilename = "../Inputs/TJ.json";
  JSONData     = ReadJSONFile (JSONFilename);
  VIZ          = JSONData["VIZ"].get<int> ();

  // ------------
  // Sanity check
  // ------------
  if (qc < 0.)
    {
      printf ("Equilibrium:: Error - qc cannot be negative\n");
      exit (1);
    }
  if (qa < qc)
    {
      printf ("Equilibrium:: Error - qa cannot be less than qc\n");
      exit (1);
    }
  if (pc < 0.)
    {
      printf ("Equilibrium:: Error - pc cannot be negative\n");
      exit (1);
    }
  if (mu < 1.)
    {
      printf ("Equilibrium:: Error - mu cannot be less than unity\n");
      exit (1);
    }
  if (epsa <= 0.)
    {
      printf ("Equilibrium:: Error - epsa must be positive\n");
      exit (1);
    }
  if (eps <= 0.)
    {
      printf ("Equilibrium:: Error - eps must be positive\n");
      exit (1);
    }
  if (Ns < 1)
    {
      printf ("Equilibrium:: Error - Ns cannot be less than unity\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("Equilibrium:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (Nf < 2)
    {
      printf ("Equilibrium:: Error - Nf cannot be less than two\n");
      exit (1);
    }
  if (Nw < 2)
    {
      printf ("Equilibrium:: Error - Nw cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("Equilibrium:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("Equilibrium:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("Equilibrium:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("Equilibrium:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("Equilibrium:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (Hna.size() != Vna.size())
    {
      printf ("Equilibrium:: Error - Hna and Van arrays must be the same size\n");
      exit (1);
    }
  if (tilt < -180. || tilt > 180.)
    {
       printf ("Equilibrium:: Error -  tilt must lie in range -180 to +180\n");
      exit (1);
    }
}

// ##########
// Destructor
// ##########
Equilibrium::~Equilibrium ()
{
}

// ##################################################################################
// Function to find nu value that gives target edge safety-factor read from JSON file
// ##################################################################################
void Equilibrium::Setnu ()
{
  LightEquilibrium lightequilibrium (qc, epsa, pc, Hna, Vna);
      
  lightequilibrium.GetNu (qa, nu);
}

// #####################################
// Function to solve equilibrium problem
// #####################################
void Equilibrium::Solve ()
{
  // .............................
  // Output calculation parameters
  // .............................
  printf ("\n");
  printf ("Class EQUILIBRIUM::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Calculation parameters:\n");
  printf ("qc  = %10.3e qa   = %10.3e epsa = %10.3e pc = %10.3e mu = %10.3e Ns = %3d Nr = %3d Nf = %3d Nw = %3d\n",
	  qc, qa, epsa, pc, mu, Ns, Nr, Nf, Nw);
  printf ("SRC =  %1d         tilt = %10.3e\n",
	  SRC, tilt);

  // ....................
  // Read in profile data
  // ....................
  if (SRC)
    {
      FILE* file = OpenFiler ("../Inputs/Profile.txt");
      if (fscanf (file, "%d", &NPTS) != 1)
	{
	  printf ("Equilibrium::Solve: Error reading NPTS from Profile.txt\n");
	  exit (1);
	}
      if (NPTS < 1)
	{
	  printf ("Equilibrium::Solve: NPTS must be greater than unity\n");
	}

      rin  = new double[NPTS];
      qin  = new double[NPTS];
      f1in = new double[NPTS];
      p2in = new double[NPTS];

      for (int i = 0; i < NPTS; i++)
	if (fscanf (file, "%lf %lf %lf", &rin[i], &qin[i], &p2in[i]) != 3)
	  {
	    printf ("Equilibrium::Solve: Error reading data from Profile.txt\n");
	    exit (1);
	  }
      
      fclose (file);

      qc = qin[0];
      qa = qin[NPTS-1];
      nu = qa /qc;
      
      double p0 = p2in[0];
      for (int i = 0; i < NPTS; i++)
	{
	  f1in[i] = rin[i]*rin[i] /qin[i];
	  p2in[i] = pc * p2in[i] /p0;
	}

      f1inspline = gsl_spline_alloc (gsl_interp_cspline, NPTS);
      p2inspline = gsl_spline_alloc (gsl_interp_cspline, NPTS);

      f1inacc = gsl_interp_accel_alloc ();
      p2inacc = gsl_interp_accel_alloc ();

      gsl_spline_init (f1inspline, rin, f1in, NPTS);
      gsl_spline_init (p2inspline, rin, p2in, NPTS);
    }

  // ...............
  // Allocate memory
  // ...............
  rr    = new double[Nr+1];
  p2    = new double[Nr+1];
  f1    = new double[Nr+1];
  f3    = new double[Nr+1];
  g2    = new double[Nr+1];
  q0    = new double[Nr+1];
  q2    = new double[Nr+1];
  It    = new double[Nr+1];
  Ip    = new double[Nr+1];
  Jt    = new double[Nr+1];
  Jp    = new double[Nr+1];
  pp    = new double[Nr+1];
  ppp   = new double[Nr+1];
  qq    = new double[Nr+1];
  qqq   = new double[Nr+1];
  s     = new double[Nr+1];
  s2    = new double[Nr+1];
  s0    = new double[Nr+1];
  S1    = new double[Nr+1];
  S2    = new double[Nr+1];
  S3    = new double[Nr+1];
  S4    = new double[Nr+1];
  S5    = new double[Nr+1];
  P1    = new double[Nr+1];
  P1a   = new double[Nr+1];
  P2    = new double[Nr+1];
  P2a   = new double[Nr+1];
  P3    = new double[Nr+1];
  P3a   = new double[Nr+1];
  ff    = new double[Nr+1];
  ggr2  = new double[Nr+1];
  RR2   = new double[Nr+1];
  IR2   = new double[Nr+1];
  Psi   = new double[Nr+1];
  PsiN  = new double[Nr+1];
  Tf    = new double[Nr+1];
  mu0P  = new double[Nr+1];
  DI    = new double[Nr+1];
  DR    = new double[Nr+1];
  Te    = new double[Nr+1];
  ne    = new double[Nr+1];
  Tep   = new double[Nr+1];
  nep   = new double[Nr+1];
  Lfunc = new double[Nr+1];

  req    = new double[2*Nf];
  weq    = new double[2*Nf];
  teq    = new double[2*Nf];
  Req    = new double[2*Nf];
  Zeq    = new double[2*Nf];
  BReq   = new double[2*Nf];
  neeq   = new double[2*Nf];
  Teeq   = new double[2*Nf];
  dRdreq = new double[2*Nf];
  dRdteq = new double[2*Nf];
  dZdreq = new double[2*Nf];
  dZdteq = new double[2*Nf];
  Weeq   = new double[2*Nf];
  wLeq   = new double[2*Nf];
  wUeq   = new double[2*Nf];
  wUHeq  = new double[2*Nf];

  HHfunc.resize (Ns+1, Nr+1);
  VVfunc.resize (Ns+1, Nr+1);
  HPfunc.resize (Ns+1, Nr+1);
  VPfunc.resize (Ns+1, Nr+1);

  Itspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Ipspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  g2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  fspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  q2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  gr2spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  R2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1); 
  I2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  sspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  qspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Lspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  wspline   = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Rspline   = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Zspline   = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);

  Itacc  = gsl_interp_accel_alloc ();
  Ipacc  = gsl_interp_accel_alloc ();
  g2acc  = gsl_interp_accel_alloc ();
  facc   = gsl_interp_accel_alloc ();
  q2acc  = gsl_interp_accel_alloc ();
  gr2acc = gsl_interp_accel_alloc ();
  R2acc  = gsl_interp_accel_alloc ();
  I2acc  = gsl_interp_accel_alloc ();
  sacc   = gsl_interp_accel_alloc ();
  qacc   = gsl_interp_accel_alloc ();
  Lacc   = gsl_interp_accel_alloc ();
  wacc   = gsl_interp_accel_alloc ();
  Racc   = gsl_interp_accel_alloc ();
  Zacc   = gsl_interp_accel_alloc ();

  HHspline = new gsl_spline*[Ns+1];
  VVspline = new gsl_spline*[Ns+1];
  HPspline = new gsl_spline*[Ns+1];
  VPspline = new gsl_spline*[Ns+1];

  HHacc = new gsl_interp_accel*[Ns+1];
  VVacc = new gsl_interp_accel*[Ns+1];
  HPacc = new gsl_interp_accel*[Ns+1];
  VPacc = new gsl_interp_accel*[Ns+1];

  for (int n = 0; n <= Ns; n++)
    {
      HHspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VVspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      HPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

      HHacc[n] = gsl_interp_accel_alloc ();
      VVacc[n] = gsl_interp_accel_alloc ();
      HPacc[n] = gsl_interp_accel_alloc ();
      VPacc[n] = gsl_interp_accel_alloc ();
     }

  RR    .resize (Nf, Nw+1);
  ZZ    .resize (Nf, Nw+1);
  dRdr  .resize (Nf, Nw+1);
  dRdt  .resize (Nf, Nw+1);
  dZdr  .resize (Nf, Nw+1);
  dZdt  .resize (Nf, Nw+1);
  Jac   .resize (Nf, Nw+1);
  Jax   .resize (Nf, Nw+1);
  RRw   .resize (Nf, Nw+1);
  ZZw   .resize (Nf, Nw+1);
  rvals .resize (Nf, Nw+1);
  thvals.resize (Nf, Nw+1);
  wvals .resize (Nf, Nw+1);

  Rbound   = new double[Nw+1];
  Zbound   = new double[Nw+1];
  tbound   = new double[Nw+1];
  wbound0  = new double[Nw+1];
  tbound0  = new double[Nw+1];
  wbound   = new double[Nw+1];
  R2b      = new double[Nw+1];
  grr2b    = new double[Nw+1];
  dRdtheta = new double[Nw+1];
  dZdtheta = new double[Nw+1];

  Rwall    = new double[Nw+1];
  Zwall    = new double[Nw+1];
  wwall    = new double[Nw+1];
  twall    = new double[Nw+1];
  R2w      = new double[Nw+1];
  grr2w    = new double[Nw+1];
  dRdthetw = new double[Nw+1];
  dZdthetw = new double[Nw+1];

  // ..................
  // Set up radial grid
  // ..................
  for (int i = 0; i <= Nr; i++)
    {
      rr [i] = double (i) /double (Nr);
      p2 [i] = Getp2  (rr[i]);
      pp [i] = Getp2p (rr[i]);
      ppp[i] = Getp2pp(rr[i]);
      f1 [i] = Getf1  (rr[i]);
    }

  // ....................................
  // Integrate shaping function equations
  // ....................................
  double  r, h, t_err;
  int     rept;
  int     neqns = 1 + 2 + 4 * (Ns - 1);
  double* y     = new double[neqns];
  double* err   = new double[neqns];
  rhs_chooser   = 0;

  h            = h0;
  count        = 0;
  r            = rr[0];
  g2[0]        = 0.;
  HHfunc(1, 0) = 0.;
  HPfunc(1, 0) = 0.;
  for (int n = 2; n <= Ns; n++)
    {
      if (n == 2)
	{
	  HHfunc(n, 0) = 0.;
	  HPfunc(n, 0) = 1.;
	  VVfunc(n, 0) = 0.;
	  VPfunc(n, 0) = 1.;
	}
      else
	{
	  HHfunc(n, 0) = 0.;
	  HPfunc(n, 0) = 0.;
	  VVfunc(n, 0) = 0.;
	  VPfunc(n, 0) = 0.;
	}
    }

  double f1c   = 1./qc;
  double p2ppc = - 2.*mu*pc;

  r     = eps;
  y[0]  = - (f1c*f1c + p2ppc/2.)  *r*r;
  y[1]  = (2.*p2ppc/f1c/f1c - 1.) *r*r/8.;
  y[2]  = (2.*p2ppc/f1c/f1c - 1.) *r  /4.;
  int j = 3;
  for (int n = 2; n <= Ns; n++)
    {
      y[j] =                pow (r, double (n-1)); j++;
      y[j] = double (n-1) * pow (r, double (n-2)); j++;
      y[j] =                pow (r, double (n-1)); j++;
      y[j] = double (n-1) * pow (r, double (n-2)); j++;
    }

  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (neqns, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      CashKarp45Fixed (neqns, r, y, err, rr[i] - r);

      g2[i]        = y[0];
      HHfunc(1, i) = y[1];
      HPfunc(1, i) = y[2];
      int j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  HHfunc(n, i) = y[j]; j++;
	  HPfunc(n, i) = y[j]; j++;
	  VVfunc(n, i) = y[j]; j++;
	  VPfunc(n, i) = y[j]; j++;
	}
    }

  delete[] y; delete[] err;

  // .....................
  // Set edge shaping data
  // .....................
  int     nshape = Hna.size();
  double* hna    = new double[nshape+2];
  double* vna    = new double[nshape+2];

  for (int i = 0; i < nshape; i++)
    {
      hna[i+2] = Hna[i];
      vna[i+2] = Vna[i];
     }

  if (nshape > Ns)
    nshape = Ns;
  
  double f1a = f1[Nr];
  double H1a = HPfunc(1, Nr);
  
  int    nn   = 1;
  double zero = 0.;
  printf ("n = %3d:  Hna = %10.3e  Vna = %10.3e\n", nn, HHfunc(1, Nr), zero);

  // .........................
  // Rescale shaping functions
  // .........................
  for (int n = 2; n <= Ns; n++)
    {
      double Hnfc, Vnfc, Hnfca, Vnfca;
      if (n <= nshape+1)
	{
	  Hnfc  = hna[n];
	  Vnfc  = vna[n];
	  Hnfca = HHfunc(n, Nr);
	  Vnfca = VVfunc(n, Nr);
	}
      else
	{
	  Hnfc  = 0.;
	  Vnfc  = 0.;
	  Hnfca = 1.;
	  Vnfca = 1.;
	}
      
      for (int i = 0; i <= Nr; i++)
	{
	  HHfunc(n, i) = HHfunc(n, i) * Hnfc /Hnfca; 
	  HPfunc(n, i) = HPfunc(n, i) * Hnfc /Hnfca; 
	  VVfunc(n, i) = VVfunc(n, i) * Vnfc /Vnfca; 
	  VPfunc(n, i) = VPfunc(n, i) * Vnfc /Vnfca; 
	}

      double Hnam = (HPfunc(n, Nr) + double (n - 1) * HHfunc(n, Nr)) /double (2*n);
      double Vnam = (VPfunc(n, Nr) + double (n - 1) * VVfunc(n, Nr)) /double (2*n);

      if (fabs (Hnfc) > 1.e-15 || fabs (Vnfc) > 1.e-15)
	printf ("n = %3d:  Hna = %10.3e  Vna = %10.3e\n", n, Hnfc, Vnfc);
    }
  
  delete[] hna; delete[] vna;

  // .............................
  // Interpolate shaping functions
  // .............................
  gsl_spline_init (g2spline, rr, g2, Nr+1);

  double* data = new double[Nr+1];

  for (int i = 0; i <= Nr; i++)
    data[i] = HHfunc(1, i);
  gsl_spline_init (HHspline[1], rr, data, Nr+1);

  for (int i = 0; i <= Nr; i++)
    data[i] = HPfunc(1, i);
  gsl_spline_init (HPspline[1], rr, data, Nr+1);

  for (int n = 2; n <= Ns; n++)
    {
      for (int i = 0; i <= Nr; i++)
	data[i] = HHfunc(n, i);
      gsl_spline_init (HHspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = HPfunc(n, i);
      gsl_spline_init (HPspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VVfunc(n, i);
      gsl_spline_init (VVspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VPfunc(n, i);
      gsl_spline_init (VPspline[n], rr, data, Nr+1);
    }

  delete[] data;
  
  // .................................
  // Interpolate relabelling parameter
  // .................................
  Lfunc[0] = 0.;
  for (int i = 1; i <= Nr; i++)
    {
      double rf = double (i) /double (Nr);

      Lfunc[i] = GetL (rf, 1);
    }

  gsl_spline_init (Lspline, rr, Lfunc, Nr+1);

  // .....................
  // Integrate f3 equation
  // .....................
  printf ("Calculating f3 equation:\n");
  double* y1    = new double[2];
  double* dy1dr = new double[2];
  double* err1  = new double[2];
  rhs_chooser   = 1;

  f1c        = 1. /qc;
  double H2c = HPfunc(2, 0);
  double V2c = VPfunc(2, 0);

  h       = h0;
  count   = 0;
  f3[0]   = 0.;
  Psi[0]  = 0.;
  q0[0]   = qc;
  qq[0]   = 0.;
  ff[0]   = 0.;
  ggr2[0] = 1. + epsa*epsa * (H2c*H2c + V2c*V2c);
  RR2[0]  = 1.;
  IR2[0]  = 1.;
  It[0]   = 0.;
  Ip[0]   = 0.;

  r     = eps;
  y1[0] = - f1c * (H2c*H2c + V2c*V2c)   * r*r;
  y1[1] = 0.5 * (f1c * eps*eps * y1[0]) * r*r;
  
  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      CashKarp45Fixed (2, r, y1, err1, rr[i] - r);
      CashKarp45Rhs (rr[i], y1, dy1dr);
      
      f3[i]  = y1[0];
      Psi[i] = y1[1];
      q0[i]  = rr[i]*rr[i] /f1[i];
      q2[i]  = rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) * exp(- epsa*epsa * f3[i]/f1[i]) /f1[i];

      double ff1 = f1[i];
      double ff3 = f3[i];
      double f1p = Getf1p (r);
      double p2p = pp[i];
      double gg2 = g2[i];
      double g2p = - p2p - ff1*f1p/r/r;
      double f3p = dy1dr[0];

      qq[i] =     2.*rr[i]*rr[i] * (1. + epsa*epsa*g2[i])                                           * exp(- epsa*epsa * f3[i]/f1[i]) /f1[i]
	     + rr[i]*rr[i]*rr[i] * (     epsa*epsa*g2p  )                                           * exp(- epsa*epsa * f3[i]/f1[i]) /f1[i]
	     - rr[i]*rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) * epsa*epsa * (f3p/ff1 - f1p*ff3/ff1/ff1) * exp(- epsa*epsa * f3[i]/f1[i]) /f1[i]
	     - rr[i]*rr[i]*rr[i] * (1. + epsa*epsa*g2[i])                                           * exp(- epsa*epsa * f3[i]/f1[i]) * f1p/ff1/ff1;
  
      double gr2 = 3. *rr[i]*rr[i]/4. -    HHfunc(1, i)                  + HPfunc(1, i) * HPfunc(1, i) /2.;
      double ir2 = 13.*rr[i]*rr[i]/4. - 3.*HHfunc(1, i) + r*HPfunc(1, i) + HPfunc(1, i) * HPfunc(1, i) /2.;
      for (int n = 2; n <= Ns; n++)
	{
	  gr2 += (HPfunc(n, i) * HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i) * HHfunc(n, i) /r/r) /2.;
	  gr2 += (VPfunc(n, i) * VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i) * VVfunc(n, i) /r/r) /2.;
	  ir2 += (HPfunc(n, i) * HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i) * HHfunc(n, i) /r/r) /2.;
	  ir2 += (VPfunc(n, i) * VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i) * VVfunc(n, i) /r/r) /2.;
	}

      double R2 = rr[i]*rr[i]/2. - rr[i]*HPfunc(1, i) - 2.*HHfunc(1, i);

      ff  [i] = ff1 + epsa*epsa * ff3;
      ggr2[i] = 1.  + epsa*epsa * gr2;
      RR2 [i] = 1.  - epsa*epsa * R2;
      IR2 [i] = 1.  + epsa*epsa * ir2;
      It  [i] =   2.*M_PI * (f1[i] + epsa*epsa * f3[i] + epsa*epsa * f1[i] * gr2);
      Ip  [i] = - 2.*M_PI * g2[i];
    }
  q2[0] = qc * (1. + epsa*epsa * (H2c*H2c + V2c*V2c));

  gshape = It[Nr] /2./M_PI /f1[Nr];

  for (int i = 0; i <= Nr; i++)
    {
      Tf[i]   = B0*R0 * (1. + epsa*epsa*g2[i]);
      mu0P[i] = B0*B0 * epsa*epsa * p2[i];
      PsiN[i] = Psi[i] /Psi[Nr];
    }

  delete[] y1; delete[] dy1dr; delete[] err1;

  // ........................
  // Calculate magnetic shear
  // ........................
  for (int i = 0; i <= Nr; i++)
    {
      s[i] = qq[i] /q2[i];
    }

  for (int i = 1; i <= Nr; i++)
    {
      s0[i] = 2. - rr[i] * Getf1p (rr[i]) /Getf1 (rr[i]);
    }
  s0[0] = 0.;

  // ....................................
  // Calculate target edge magnetic shear
  // ....................................
  double sum  = 1.5 - 2. * HPfunc(1, Nr) + HPfunc(1, Nr) * HPfunc(1, Nr);
  for (int n = 2; n <= Ns; n++)
    {
      sum +=
	+ HPfunc(n, Nr) * HPfunc(n, Nr) + 2. * double (n*n - 1) * HPfunc(n, Nr) * HHfunc(n, Nr) - double (n*n - 1) * HHfunc(n, Nr) * HHfunc(n, Nr)
	+ VPfunc(n, Nr) * VPfunc(n, Nr) + 2. * double (n*n - 1) * VPfunc(n, Nr) * VVfunc(n, Nr) - double (n*n - 1) * VVfunc(n, Nr) * VVfunc(n, Nr);
    }
  double sa = 2. + epsa*epsa * sum;
   
  // ........................
  // Calculate s2 = r^2 q''/q
  // ........................
  double hh = 1. /double(Nr);
  for (int i = 1; i < Nr; i++)
    qqq[i] = rr[i] * (qq[i+1] - qq[i-1]) /2./hh;
  qqq[0]  = 0.;
  qqq[Nr] = qqq[Nr-1] + (qqq[Nr-1] - qqq[Nr-2]);

  for (int i = 0; i <= Nr; i++)
    {
      s2[i] = (qqq[i] - qq[i]) /q2[i];
    }

  // ...............................
  // Interpolate equilibrium splines
  // ...............................
  gsl_spline_init (fspline,   rr, ff,   Nr+1);
  gsl_spline_init (q2spline,  rr, q2,   Nr+1);
  gsl_spline_init (gr2spline, rr, ggr2, Nr+1);
  gsl_spline_init (R2spline,  rr, RR2,  Nr+1);
  gsl_spline_init (I2spline,  rr, IR2,  Nr+1);
  gsl_spline_init (sspline,   rr, s,    Nr+1);
  gsl_spline_init (qspline,   rr, q2,   Nr+1);
  gsl_spline_init (Itspline,  rr, It,   Nr+1);
  gsl_spline_init (Ipspline,  rr, Ip,   Nr+1);

  for (int i = 0; i <= Nr; i++)
    {
      Jt [i] = gsl_spline_eval_deriv (Itspline, rr[i], Itacc);
      Jp [i] = gsl_spline_eval_deriv (Ipspline, rr[i], Ipacc);
      DI [i] = GetDI  (rr[i]);
      DR [i] = GetDR  (rr[i]);
      ne [i] = Getne  (rr[i]);
      Te [i] = GetTe  (rr[i]);
      nep[i] = Getnep (rr[i]);
      Tep[i] = GetTep (rr[i]);
    }
  DI[0] = DI[1];
  DR[0] = DR[1];
  ne [Nr] = 2.*ne [Nr-1] - ne [Nr-2];
  Te [Nr] = 2.*Te [Nr-1] - Te [Nr-2];
  nep[Nr] = 2.*nep[Nr-1] - nep[Nr-2];
  Tep[Nr] = 2.*Tep[Nr-1] - Tep[Nr-2];

  // ............................
  // Calculate li and beta values
  // ............................
  double* y2   = new double[4];
  double* err2 = new double[4];
  rhs_chooser  = 2;

  r     = eps;
  h     = h0;
  count = 0;
  y2[0] = 0.;
  y2[1] = 0.; 
  y2[2] = 0.;
  y2[3] = 0.;

  do
    {
      CashKarp45Adaptive (4, r, y2, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < 1. - h);
  CashKarp45Fixed (4, r, y2, err2, 1. - r);

  amean = (GetR (1., M_PI, 1) - GetR (1., 0., 1)) /2.;
  li    = 2. * y2[0] /ff[Nr]/ff[Nr] /ggr2[Nr] /ggr2[Nr];
  betat = 2. * epsa*epsa * y2[1] /y2[3];
  betap = 2. * y2[1] /y2[3] /ff[Nr]/ff[Nr] /ggr2[Nr] /ggr2[Nr];
  betaN = 20. * betat * (amean /epsa) /epsa /ff[Nr] /ggr2[Nr];

  delete[] y2; delete[] err2;

  printf ("qc = %10.3e q0a   = %10.3e q2a   = %10.3e Ip    = %10.3e It = %10.3e gshape = %10.3e\n",
	  q2[0], q0[Nr], q2[Nr], Ip[Nr], It[Nr], gshape);
  printf ("li = %10.3e betat = %10.3e betap = %10.3e betaN = %10.3e\n",
  	  li, betat, betap, betaN);

  // ...........................
  // Calculate shaping functions
  // ...........................
  for (int i = 1; i <= Nr; i++)
    {
      double sum1 = 1.5 * HPfunc(1, i) * HPfunc(1, i);
      double sum2 = 1.5 * rr[i]*rr[i] - 2. * rr[i] * HPfunc(1, i) + HPfunc(1, i)*HPfunc(1, i);
      double sum3 = - 0.75 * rr[i]*rr[i] + rr[i]*rr[i] /q2[i]/q2[i] + HHfunc(1, i) + 1.5 * HPfunc(1, i)*HPfunc(1, i);
      double sum4 = 0.75 * rr[i]*rr[i] - HHfunc(1, i) + HPfunc(1, i) * HPfunc(1, i) /2.;
      double sum5 = 1.75 * rr[i]*rr[i] - HHfunc(1, i) - 3. * rr[i] * HPfunc(1, i) + 1.5 * HPfunc(1, i) * HPfunc(1, i);
 
      for (int n = 2; n <= Ns; n++)
	{
	  sum1 += (                    3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		       - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]) /2.;
	  sum2 +=                            HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i)
	          + 2. * double (n*n - 1) * (HPfunc(n, i) * HHfunc(n, i) + VPfunc(n, i) * VVfunc(n, i)) /rr[i]
	          -      double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i];
	  sum3 += (                    3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		       - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]) /2.;
	  sum4 +=                           (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i)) /2.
	               + double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i] /2.;
	  sum5 += (                    3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		       - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]) /2.;
	}

      S1[i] = sum1;
      S2[i] = sum2;
      S3[i] = sum3;
      S4[i] = sum4;
      S5[i] = sum5;
    }
  S1[0] = 0.;
  S2[0] = S2[1];
  S3[0] = 0.;
  S4[0] = S4[1];
  S5[0] = 0.;

  // ...........................
  // Calculate profile functions
  // ...........................
  for (int i = 0; i <= Nr; i++)
    {
      P1 [i] = (2. - s[i]) /q2[i];
      P1a[i] = (2. - s0[i]) /q0[i];
      P2 [i] = (- 3. * s[i] + 2. * s[i]*s[i] - s2[i]) /q2[i];
      P3a[i] = - (2. - s[i]) * S3[i] /q2[i] + S2[i] /q2[i];
  }

  for (int i = 1; i < Nr; i++)
    {
      P2a[i] = rr[i] * (P1a[i+1] - P1a[i-1]) /2./hh;
      P3[i]  = 2. * rr[i] * pp[i] * (2. - s[i]) - q2[i] * rr[i] * (P3a[i+1] - P3a[i-1]) /2./hh;
    }
  P2a[0]  = P2a[1]    - (P2a[2]    - P2a[1]);
  P2a[Nr] = P2a[Nr-1] + (P2a[Nr-1] - P2a[Nr-2]);

  P3[0]  = P3[1]    - (P3[2]    - P3[1]);
  P3[Nr] = P3[Nr-1] + (P3[Nr-1] - P3[Nr-2]);

  if (VIZ)
    {
      // ...........................................................
      // Calculate magnetic flux-surfaces for visualization purposes
      // ...........................................................
      printf ("Calculating magnetic flux-surfaces:\n");
      
      for (int i = 1; i <= Nf; i++)
	{
	  double rf = double (i) /double (Nf);
	  
	  for (int j = 0; j <= Nw; j++)
	    {
	      double t = double (j) * 2.*M_PI /double (Nw);
	      
	      double w, wold = t, dtdw, dtdwold = 1.;
	      for (int ii = 0; ii < 5; ii++)
		{
		  w    = t - Gettheta (rf, wold, 1);
		  wold = w;

		  dtdw    = 1. + Getdthetadomega (rf, w, 1);
		  dtdwold = dtdw;
		}
	      
	      double R  = GetR    (rf, w, 1);
	      double Z  = GetZ    (rf, w, 1);
	      double Rr = GetdRdr (rf, w, 1);
	      double Rw = GetdRdw (rf, w, 1);
	      double Zr = GetdZdr (rf, w, 1);
	      double Zw = GetdZdw (rf, w, 1);
	      
	      RR    (i-1, j) = R;
	      ZZ    (i-1, j) = Z;
	      dRdr  (i-1, j) = Rr;
	      dRdt  (i-1, j) = Rw /dtdw /rf;
	      dZdr  (i-1, j) = Zr;
	      dZdt  (i-1, j) = Zw /dtdw /rf;
	      Jac   (i-1, j) = (Rw * Zr - Rr * Zw) /dtdw /rf;
	      Jax   (i-1, j) = R;
	      rvals (i-1, j) = rf;
	      thvals(i-1, j) = t;
	      wvals (i-1, j) = w;
	    }
	}
      
      for (int i = 1; i <= Nf; i++)
	{
	  double rf = double (i) /double (Nf);
	  
	  for (int j = 0; j <= Nw; j++)
	    {
	      double w = double (j) * 2.*M_PI /double (Nw);
	      
	      double R = GetR (rf, w, 1);
	      double Z = GetZ (rf, w, 1);
	      
	      RRw(i-1, j) = R;
	      ZZw(i-1, j) = Z;
	    }
	}

      // ...................................
      // Calculate synthetic diagnostic data
      // ...................................
      for (int i = 0; i < Nf; i++)
	{
	  double rf = 1. - double (i) /double (Nf);
	  double wt = tilt*M_PI/180.;
	  double wf = wt;

	  double Ro, Z0, dZ0dw;
	  for (int j = 0; j < 5; j++)
	    {
	      Ro    = GetR    (rf, wf, 1);
	      Z0    = GetZ    (rf, wf, 1);
	      dZ0dw = GetdZdw (rf, wf, 1);

	      wf += ((1. - Ro) * tan(wt) - Z0) /dZ0dw;
	    }

	  req[i] = rf;
	  weq[i] = wf;
	  teq[i] = wf + Gettheta (rf, wf, 1);
	  Req[i] = GetR (rf, wf, 1);
	  Zeq[i] = GetZ (rf, wf, 1);

	  if (teq[i] < 0.)
	    teq[i] += 2.*M_PI;
	  if (teq[i] > 2.*M_PI)
	    teq[i] -= 2.*M_PI;

	  double Rr   = GetdRdr (rf, wf, 1);
	  double Rw   = GetdRdw (rf, wf, 1);
	  double Zr   = GetdZdr (rf, wf, 1);
	  double Zw   = GetdZdw (rf, wf, 1);
	  double dtdw = 1. + Getdthetadomega (rf, wf, 1);
	  double g2   = gsl_spline_eval (g2spline, rf, g2acc);
	  double q    = Getq (rf);
	      
	  dRdreq[i] = Rr;
	  dRdteq[i] = Rw /dtdw /rf;
	  dZdreq[i] = Zr;
	  dZdteq[i] = Zw /dtdw /rf;

	  BReq[i] = epsa * B0 * (1. + epsa*epsa*g2) * (Rw * cos(wt) - Zw * sin(wt)) /q /Req[i]/Req[i] /dtdw;
	  neeq[i] = Getne (rf);
	  Teeq[i] = GetTe (rf);

	  double BT = B0 /Req[i];
	  double We = pc_e * BT /pc_m_e;
	  double wp = sqrt (neeq[i] * pc_e*pc_e /pc_epsilon_0/pc_m_e);

	  Weeq [i] = We;
	  wUeq [i] =   0.5 * We + sqrt (0.25*We*We + wp*wp);
	  wLeq [i] = - 0.5 * We + sqrt (0.25*We*We + wp*wp);
	  wUHeq[i] = sqrt (We*We + wp*wp);
	}
      neeq[0] = 2.*neeq[1] - neeq[2];
      Teeq[0] = 2.*Teeq[1] - Teeq[2];

      for (int i = 0; i < Nf; i++)
	{
	  double rf = double (i+1) /double (Nf);
	  double wt = M_PI + tilt*M_PI/180.;
	  double wf = wt;

	  double Ro, Z0, dZ0dw;
	  for (int j = 0; j < 5; j++)
	    {
	      Ro    = GetR    (rf, wf, 1);
	      Z0    = GetZ    (rf, wf, 1);
	      dZ0dw = GetdZdw (rf, wf, 1);
	      
	      wf += ((1. - Ro) * tan(wt) - Z0) /dZ0dw;
	    }

	  req[i+Nf] = rf;
	  weq[i+Nf] = wf;
	  teq[i+Nf] = wf + Gettheta (rf, wf, 1);
	  Req[i+Nf] = GetR (rf, wf, 1);
	  Zeq[i+Nf] = GetZ (rf, wf, 1);

	  if (teq[i] < 0.)
	    teq[i] += 2.*M_PI;
	  if (teq[i] > 2.*M_PI)
	    teq[i] -= 2.*M_PI;

	  double Rr   = GetdRdr (rf, wf, 1);
	  double Rw   = GetdRdw (rf, wf, 1);
	  double Zr   = GetdZdr (rf, wf, 1);
	  double Zw   = GetdZdw (rf, wf, 1);
	  double dtdw = 1. + Getdthetadomega (rf, wf, 1);
	  double g2   = gsl_spline_eval (g2spline, rf, g2acc);
	  double q    = Getq (rf);
	      
	  dRdreq[i+Nf] = Rr;
	  dRdteq[i+Nf] = Rw /dtdw /rf;
	  dZdreq[i+Nf] = Zr;
	  dZdteq[i+Nf] = Zw /dtdw /rf;

	  BReq[i+Nf] = epsa * B0 * (1. + epsa*epsa*g2) * (Rw * cos(wt) - Zw * sin(wt)) /q /Req[i+Nf]/Req[i+Nf] /dtdw;
	  neeq[i+Nf] = Getne (rf);
	  Teeq[i+Nf] = GetTe (rf);

	  double BT = B0 /Req[i+Nf];
	  double We = pc_e * BT /pc_m_e;
	  double wp = sqrt (neeq[i+Nf] * pc_e*pc_e /pc_epsilon_0/pc_m_e);

	  Weeq [i+Nf] = We;
	  wUeq [i+Nf] =   0.5 * We + sqrt (0.25*We*We + wp*wp);
	  wLeq [i+Nf] = - 0.5 * We + sqrt (0.25*We*We + wp*wp);
	  wUHeq[i+Nf] = sqrt (We*We + wp*wp);
	}
      neeq[2*Nf-1] = 2.*neeq[2*Nf-2] - neeq[2*Nf-3];
      Teeq[2*Nf-1] = 2.*Teeq[2*Nf-2] - Teeq[2*Nf-3];
    }
  
  // .......................
  // Calculate boundary data
  // .......................
  printf ("Calculating boundary data:\n");

  // Set up preliminary omega grid
  for (int j = 0; j <= Nw; j++)
    wbound0[j] = double (j) * 2.*M_PI /double (Nw);
  
  // Calculate preliminary theta grid
  tbound0[0] = 0.;
  
  double  w;
  double* y3   = new double[1];
  double* err3 = new double[1];
  rhs_chooser  = 3;
  
  w     = 0.;
  h     = h0;
  count = 0;
  y3[0] = 0.;
  
  for (int j = 1; j <= Nw; j++)
    {
      do
	{
	  CashKarp45Adaptive (1, w, y3, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (w < wbound0[j]);
      CashKarp45Fixed (1, w, y3, err3, wbound0[j] - w);
      
      tbound0[j] = y3[0];
    }

  double raa = tbound0[Nw] /(2.*M_PI);
  printf ("ra = %10.3e\n", raa);
  
  for (int j = 0; j < Nw; j++)
    tbound0[j] /= raa;

  tbound0[Nw] = 2.*M_PI;
  
  // Interpolate preliminary boundary data
  gsl_spline_init (wspline, tbound0, wbound0, Nw+1);
  
  // Calculate final boundary data
  for (int j = 0; j <= Nw; j++)
    tbound[j] = double (j) * 2.*M_PI /double (Nw);
  
  double ra = 1.;
  for (int j = 0; j <= Nw; j++)
    {
      double t = tbound[j];
      double w = gsl_spline_eval (wspline, t, wacc);
      
      wbound[j] = w;
      Rbound[j] = GetR    (ra, w, 1);
      Zbound[j] = GetZ    (ra, w, 1);
      R2b   [j] = GetR2   (ra, t);
      grr2b [j] = Getgrr2 (ra, t);
    }
  
  gsl_spline_init (Rspline, tbound, Rbound, Nw+1);
  gsl_spline_init (Zspline, tbound, Zbound, Nw+1);
  
  for (int j = 0; j <= Nw; j++)
    {
      dRdtheta[j] = gsl_spline_eval_deriv (Rspline, tbound[j], Racc);
      dZdtheta[j] = gsl_spline_eval_deriv (Zspline, tbound[j], Zacc);
    }

  // ....................
  // Read wall parameters
  // ....................
  string JSONFilename = "../Inputs/Equilibrium.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  bw = JSONData["bw"].get<double> ();

  // ...................
  // Calculate wall data
  // ...................
  H1b  = GetHHvac (1, bw);
  H1pb = GetHPvac (1, bw);

  printf ("Calculating wall data: bw = %10.3e H1 = %10.3e H1p = %10.3e\n", bw, H1b, H1pb);

  // Calculate preliminary theta grid
  tbound0[0] = 0.;

  rhs_chooser  = 4;
  
  w     = 0.;
  h     = h0;
  count = 0;
  y3[0] = 0.;
  
  for (int j = 1; j <= Nw; j++)
    {
      do
	{
	  CashKarp45Adaptive (1, w, y3, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (w < wbound0[j]);
      CashKarp45Fixed (1, w, y3, err3, wbound0[j] - w);
      
      tbound0[j] = y3[0];
    }

  double raw = tbound0[Nw] /(2.*M_PI);
  printf ("rw = %10.3e\n", raw);
  
  for (int j = 0; j < Nw; j++)
    tbound0[j] /= raw;

  tbound0[Nw] = 2.*M_PI;
  
  // Interpolate preliminary boundary data
  gsl_spline_init (wspline, tbound0, wbound0, Nw+1);
  
  // Calculate final wall data
  for (int j = 0; j <= Nw; j++)
    twall[j] = double (j) * 2.*M_PI /double (Nw);
  
  double rb = bw;
  for (int j = 0; j <= Nw; j++)
    {
      double t = twall[j];
      double w = gsl_spline_eval (wspline, t, wacc);
      
      wwall [j] = w;
      Rwall [j] = GetR    (rb, w, 1);
      Zwall [j] = GetZ    (rb, w, 1);
      R2w   [j] = GetR2   (rb, t);
      grr2w [j] = Getgrr2 (rb, t);
    }
  
  gsl_spline_init (Rspline, twall, Rwall, Nw+1);
  gsl_spline_init (Zspline, twall, Zwall, Nw+1);
  
  for (int j = 0; j <= Nw; j++)
    {
      dRdthetw[j] = gsl_spline_eval_deriv (Rspline, twall[j], Racc);
      dZdthetw[j] = gsl_spline_eval_deriv (Zspline, twall[j], Zacc);
    }

  delete[] y3; delete[] err3;
  
  // .......................
  // Calculate RMP coil data
  // .......................
  vector<double> wcoil, icoil;
  for (const auto& number : JSONData["wcoil"])
    {
      wcoil.push_back (number.get<double> ());
    }
  for (const auto& number : JSONData["Icoil"])
    {
      icoil.push_back (number.get<double> ());
    }

  if (wcoil.size() != icoil.size())
    {
      printf ("Equilibrium:: Error reading etacoil and Icoil arrays must be same size\n");
      exit (1);
    }
  ncoil = wcoil.size();

  Rcoil = new double[ncoil];
  Zcoil = new double[ncoil];
  Icoil = new double[ncoil];
  
  for (int i = 0; i < ncoil; i++)
    {
      Rcoil[i] = 1. + bw * (GetR (ra, wcoil[i]*M_PI, 1) - 1.);
      Zcoil[i] = bw * GetZ (ra, wcoil[i]*M_PI, 1);
      Icoil[i] = icoil[i];
     }
  
  printf ("RMP coil data:\n");
  for (int i = 0; i < ncoil; i++)
    printf ("w_coil/M_PI = %10.3e Rcoil = %10.3e Zcoil = %10.3e Icoil = %10.3e\n",
	    wcoil[i], Rcoil[i], Zcoil[i], Icoil[i]);

  // ......................................
  // Output equilibrium data to netcdf file
  // ......................................
  WriteNetcdf (sa);

  // ...................
  // Calculate EFIT data
  // ...................
  CalculateEFIT ();

  // ...............
  // Write EFIT file
  // ...............
  int serr = system ("../bin/write_efit");
  if (serr != 0)
    printf ("TJ:: Warning: Error runing WriteEFIT\n");
  
  // ........
  // Clean up
  // ........
  printf ("Cleaning up:\n");
  
  delete[] rr;  delete[] p2;  delete[] f1;  delete[] f3;  delete[] g2;
  delete[] q0;  delete[] q2;  delete[] It;  delete[] Ip;  delete[] Jt;
  delete[] Jp;  delete[] pp;  delete[] ppp; delete[] qq;  delete[] qqq;
  delete[] s;   delete[] s2;  delete[] S1;  delete[] S2;  delete[] P1;
  delete[] P2;  delete[] P3;  delete[] P3a; delete[] ff;  delete[] ggr2;
  delete[] RR2; delete[] IR2; delete[] S3;  delete[] P1a; delete[] P2a;
  delete[] s0;  delete[] S4;  delete[] S5;

  delete[] req;    delete[] weq;    delete[] teq;    delete[] Req; 
  delete[] Zeq;    delete[] BReq;   delete[] neeq;   delete[] Teeq;
  delete[] dRdreq; delete[] dRdteq; delete[] dZdreq; delete[] dZdteq;
  delete[] Weeq;   delete[] wUeq;   delete[] wLeq;   delete[] wUHeq;
 
  gsl_spline_free (Itspline);
  gsl_spline_free (Ipspline);
  gsl_spline_free (g2spline);
  gsl_spline_free (fspline);
  gsl_spline_free (q2spline);
  gsl_spline_free (gr2spline);
  gsl_spline_free (R2spline);
  gsl_spline_free (I2spline);
  gsl_spline_free (sspline);
  gsl_spline_free (qspline);
  gsl_spline_free (Lspline); 
  gsl_spline_free (wspline); 
  gsl_spline_free (Rspline); 
  gsl_spline_free (Zspline);

  gsl_interp_accel_free (g2acc);
  gsl_interp_accel_free (Itacc);
  gsl_interp_accel_free (Ipacc);
  gsl_interp_accel_free (facc);
  gsl_interp_accel_free (q2acc);
  gsl_interp_accel_free (gr2acc);
  gsl_interp_accel_free (R2acc);
  gsl_interp_accel_free (I2acc);
  gsl_interp_accel_free (sacc);
  gsl_interp_accel_free (qacc);
  gsl_interp_accel_free (Lacc);
  gsl_interp_accel_free (wacc);
  gsl_interp_accel_free (Racc);
  gsl_interp_accel_free (Zacc);

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

  delete[] Rbound; delete[] Zbound; delete[] tbound; delete[] wbound0;  delete[] tbound0;
  delete[] wbound; delete[] R2b;    delete[] grr2b;  delete[] dRdtheta; delete[] dZdtheta;

  delete[] Rwall; delete[] Zwall; delete[] wwall;    delete[] twall;
  delete[] R2w;   delete[] grr2w; delete[] dRdthetw; delete[] dZdthetw;

  delete[] Psi; delete[] PsiN; delete[] Tf;    delete[] mu0P;
  delete[] DI;  delete[] DR;   delete[] Lfunc; delete[] Te;
  delete[] ne;  delete[] Tep;  delete[] nep;  

  delete[] Rcoil; delete[] Zcoil; delete[] Icoil;

  if (SRC)
    {
      delete[] rin; delete[] qin; delete[] f1in; delete[] p2in;

      gsl_spline_free (f1inspline);
      gsl_spline_free (p2inspline);

      gsl_interp_accel_free (f1inacc);
      gsl_interp_accel_free (p2inacc);
    }
 } 

// ########################
// Function to return f1(r)
// ########################
double Equilibrium::Getf1 (double r)
{
  if (SRC)
    {
      return gsl_spline_eval (f1inspline, r, f1inacc);
    }
  else
    {
      if (r < 0.1)
	return r*r * (1. - (nu - 1.)*r*r/2. + (nu - 1.)*(nu - 2.)*r*r*r*r/6. - (nu - 1.)*(nu - 2.)*(nu - 3.)*r*r*r*r*r*r/24.)/qc;
      else
	return (1. - pow (1. - r*r, nu)) /nu/qc;
    }
}

// #########################
// Function to return f1'(r)
// #########################
double Equilibrium::Getf1p (double r)
{
  if (SRC)
    {
      return gsl_spline_eval_deriv (f1inspline, r, f1inacc);
    }
  else
    {
      if (r < 0.1)
	return 2.*r * (1. - (nu - 1.)*r*r + (nu - 1.)*(nu - 2.)*r*r*r*r/2. - (nu - 1.)*(nu - 2.)*(nu - 3.)*r*r*r*r*r*r/6.)/qc;
      else
	return 2. * r * pow (1. - r*r, nu-1.) /qc;
    }
}

// ########################
// Function to return p2(r)
// ########################
double Equilibrium::Getp2 (double r)
{
  if (SRC)
    {
      return gsl_spline_eval (p2inspline, r, p2inacc);
    }
  else
    {
      return pc * pow (1. - r*r, mu);
    }
}

// #########################
// Function to return p2'(r)
// #########################
double Equilibrium::Getp2p (double r)
{
  if (SRC)
    {
      return gsl_spline_eval_deriv (p2inspline, r, p2inacc);
    }
  else
    {
      return - 2. * pc * mu * r * pow (1. - r*r, mu - 1.);
    }
}

// ##########################
// Function to return p2''(r)
// ##########################
double Equilibrium::Getp2pp (double r)
{
  if (SRC)
    {
      return gsl_spline_eval_deriv2 (p2inspline, r, p2inacc);
    }
  else
    {
      return - 2. * pc * mu * pow (1. - r*r, mu - 1.)
	+ 4. * pc * mu * (mu-1.) * r*r * pow (1. - r*r, mu - 2.);
    }
}

// ####################
// Function to return q
// ####################
double Equilibrium::Getq (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (qspline, 1., qacc);
  else
    return gsl_spline_eval (qspline, r, qacc);
}

// ####################
// Function to return s
// ####################
double Equilibrium::Gets (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (sspline, 1., qacc);
  else
    return gsl_spline_eval (sspline, r, qacc);
}

// #####################
// Function to return DI
// #####################
double Equilibrium::GetDI (double r)
{
  double pp = Getp2p (r);
  double q  = Getq (r);
  double s  = Gets (r);

  return - 0.25 - epsa*epsa * 2. * r*pp * (1. - q*q) /s/s; 
}

// #####################
// Function to return DR
// #####################
double Equilibrium::GetDR (double r)
{
  double pp  = Getp2p (r);
  double H1p = GetHnp (1, r);
  double q   = Getq (r);
  double s   = Gets (r);

  return - epsa*epsa * 2. * r*pp * (1. - q*q) /s/s - epsa*epsa * 2. * pp * q*q * H1p /s; 
}

// #####################
// Function to return ne
// #####################
double Equilibrium::Getne (double r)
{
  return n0 * pow (1. - r*r, alpha) + neped;
}

// ######################
// Function to return nep
// ######################
double Equilibrium::Getnep (double r)
{
  return - 2. * n0 * alpha * r * pow (1. - r*r, alpha-1.);
}

// #####################
// Function to return Te
// #####################
double Equilibrium::GetTe (double r)
{
  double e   = 1.602176634e-19;
  double mu0 = 4.*M_PI*1.e-7;
  
  double p2 = Getp2 (r);
  double ne = Getne (r);

  if (ne <= 0.)
    return Teped;
  else
    return (B0*B0/mu0/e/2.) * epsa*epsa * p2 /ne + Teped;
}

// ######################
// Function to return Tep
// ######################
double Equilibrium::GetTep (double r)
{
  double e    = 1.602176634e-19;
  double mu0  = 4.*M_PI*1.e-7;
  
  double p2  = Getp2  (r);
  double p2p = Getp2p (r);
  double ne  = n0 * pow (1. - r*r, alpha);
  double nep = Getnep (r);

  return (B0*B0/mu0/e/2.) * epsa*epsa * (p2p /ne - p2 * nep /ne/ne);
}

// ######################
// Function to return Hnp
// ######################
double Equilibrium::GetHnp (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (HPspline[n], 1., HPacc[n]);
  else
    return gsl_spline_eval (HPspline[n], r, HPacc[n]);
}

// #######################################
// Function to return relabeling parameter
// #######################################
double Equilibrium::GetL (double r, int order)
{
  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double L = r*r*r /8. - r*hn[1] /2.;

  // ...................
  // Lowest-order result
  // ...................
  for (int n = 2; n <= Ns; n++)
    {
      L += - double (n - 1) * hn[n] * hn[n] /2./r - double (n - 1) * vn[n] * vn[n] /2./r;
    }

  // .......................
  // Higher-order correction
  // .......................
  if (order)
    {
      L += epsa * (r*r * hn[2] /4. - r*r*r * hnp[2] /4.);
      
      for (int n = 2; n <= Ns; n++)
	{
	  if (n + 1 <= Ns)
	    L +=
	      - epsa * double (n - 1) * (hnp[n] * hn[n+1] + hnp[n+1] * hn[n]) * r /2. - double (n - 1) * hn[n] * hn[n+1] /2.
	      - epsa * double (n - 1) * (vnp[n] * vn[n+1] + vnp[n+1] * vn[n]) * r /2. - double (n - 1) * vn[n] * vn[n+1] /2.;
	}
    }
  
  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return L;
}

// ####################
// Function to return R
// ####################
double Equilibrium::GetR (double r, double w, int order)
{
  double L = GetL (r, order);

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double R = 1. - epsa * r * cos (w) + epsa*epsa*epsa * L * cos (w) + epsa*epsa * hn[1];

  for (int n = 2; n <= Ns; n++)
    R += epsa*epsa * (hn[n] * cos (double (n - 1) * w) + vn[n] * sin (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return R;
}

// #######################
// Function to return dRdr
// #######################
double Equilibrium::GetdRdr (double r, double w, int order)
{
  double Lp   = GetL (r + eps, order);
  double Lm   = GetL (r - eps, order);
  double dLdr = (Lp - Lm) /2./eps;

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double dRdr = - cos (w) + epsa*epsa * dLdr * cos (w) + epsa * hnp[1];

  for (int n = 2; n <= Ns; n++)
    dRdr += epsa * (hnp[n] * cos (double (n - 1) * w) + vnp[n] * sin (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return dRdr;
}

// #######################
// Function to return dRdw
// #######################
double Equilibrium::GetdRdw (double r, double w, int order)
{
  double L = GetL (r, order);

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
 
  double dRdw = r * sin (w) - epsa*epsa * L * sin (w);

  for (int n = 2; n <= Ns; n++)
    dRdw += epsa * double (n - 1) * ( - hn[n] * sin (double (n - 1) * w) + vn[n] * cos (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return dRdw;
}

// ####################
// Function to return Z
// ####################
double Equilibrium::GetZ (double r, double w, int order)
{
  double L = GetL (r, order);

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);

  double Z = epsa * r * sin (w) - epsa*epsa*epsa * L * sin (w);

  for (int n = 2; n <= Ns; n++)
    Z += epsa*epsa * (hn[n] * sin (double (n - 1) * w) - vn[n] * cos (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return Z;
}

// #######################
// Function to return dZdr
// #######################
double Equilibrium::GetdZdr (double r, double w, int order)
{
  double Lp   = GetL (r + eps, order);
  double Lm   = GetL (r - eps, order);
  double dLdr = (Lp - Lm) /2./eps;

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double dZdr = sin (w) - epsa*epsa * dLdr * sin (w);

  for (int n = 2; n <= Ns; n++)
    dZdr += epsa * (hnp[n] * sin (double (n - 1) * w) - vnp[n] * cos (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return dZdr;
}

// #######################
// Function to return dZdw
// #######################
double Equilibrium::GetdZdw (double r, double w, int order)
{
  double L = GetL (r, order);

  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double dZdw = r * cos (w) - epsa*epsa * L * cos (w);

  for (int n = 2; n <= Ns; n++)
    dZdw += epsa * double (n - 1) * (hn[n] * cos (double (n - 1) * w) + vn[n] * sin (double (n - 1) * w));

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return dZdw;
}

// #####################
// Function to return R2
// #####################
double Equilibrium::GetR2 (double r, double t)
{
  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double R2 = 1. - 2. * epsa * r * cos (t) - epsa*epsa * (r*r/2. - r * hnp[1] - 2. * hn[1]);

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return R2;
}

// ####################################
// Function to return H(n, r) in vacuum
// ####################################
double Equilibrium::GetHHvac (int n, double r)
{
  if (n == 1)
    {
      double H1a  = HHfunc(1, Nr);
      double H1ap = HPfunc(1, Nr);

      return H1a - 0.5 * r*r * log(r) + 0.25 * (2.*H1ap + 1.) * (r*r - 1.);
    }
  else
    {
      double Hna  = HHfunc(n, Nr);
      double Hnap = HPfunc(n, Nr);

      double Hnm = (Hnap + double (n - 1) * Hna) /double (2*n);
      double Hnp = (Hnap - double (n + 1) * Hna) /double (2*n);

      return Hnm * pow (r, double (1 + n)) - Hnp * pow (r, double (1 - n));
    }
}

// #####################################
// Function to return H'(n, r) in vacuum
// #####################################
double Equilibrium::GetHPvac (int n, double r)
{
  if (n == 1)
    {
      double H1ap = HPfunc (1, Nr);

      return - r * log(r) + H1ap * r;
    }
  else
    {
      double Hna  = HHfunc(n, Nr);
      double Hnap = HPfunc(n, Nr);

      double Hnm = (Hnap + double (n - 1) * Hna) /double (2*n);
      double Hnp = (Hnap - double (n + 1) * Hna) /double (2*n);

      return double (1 + n) * Hnm * pow (r, double (n)) - double (1 - n) * Hnp * pow (r, double (-n));
    }
}

// ####################################
// Function to return V(n, r) in vacuum
// ####################################
double Equilibrium::GetVVvac (int n, double r)
{
  double Vna  = VVfunc(n, Nr);
  double Vnap = VPfunc(n, Nr);
  
  double Vnm = (Vnap + double (n - 1) * Vna) /double (2*n);
  double Vnp = (Vnap - double (n + 1) * Vna) /double (2*n);
  
  return Vnm * pow (r, double (1 + n)) - Vnp * pow (r, double (1 - n));
}

// #####################################
// Function to return V'(n, r) in vacuum
// #####################################
double Equilibrium::GetVPvac (int n, double r)
{
  double Vna  = VVfunc(n, Nr);
  double Vnap = VPfunc(n, Nr);
  
  double Vnm = (Vnap + double (n - 1) * Vna) /double (2*n);
  double Vnp = (Vnap - double (n + 1) * Vna) /double (2*n);
  
  return double (1 + n) * Vnm * pow (r, double (n)) - double (1 - n) * Vnp * pow (r, double (-n));
 }

// ################################
// Function to return PSI in vacuum
// ################################
double Equilibrium::GetPSIvac (double r)
{
  double f1a = f1[Nr];
  double f3a = f3[Nr];

  double* hnp  = new double[Ns+1];
  double* hnm  = new double[Ns+1];
  double* vnp  = new double[Ns+1];
  double* vnm  = new double[Ns+1];

  double h1a  = HHfunc(1, Nr);
  double h1ap = HPfunc(1, Nr);
  
  for (int n = 2; n <= Ns; n++)
    {
      double hna  = HHfunc(n, Nr);
      double hnap = HPfunc(n, Nr);

      hnm[n] = (hnap + double (n - 1) * hna) /double (2*n);
      hnp[n] = (hnap - double (n + 1) * hna) /double (2*n);

      double vna  = VVfunc(n, Nr);
      double vnap = VPfunc(n, Nr);
      
      vnm[n] = (vnap + double (n - 1) * vna) /double (2*n);
      vnp[n] = (vnap - double (n + 1) * vna) /double (2*n);
    }

  double sum1 =  1. - h1ap + h1ap*h1ap;
  double sum2 = - h1ap * r*r * log(r) + 0.5 * r*r * log(r) * log (r) + 0.5 * (1. + h1ap*h1ap) * (r*r - 1.);

  for (int n = 2; n <= Ns; n++)
    {
      sum1 += double (2*n * (n + 1)) * hnm[n]*hnm[n] + double (2*n * (n - 1)) * hnp[n]*hnp[n]
  	    + double (2*n * (n + 1)) * vnm[n]*vnm[n] + double (2*n * (n - 1)) * vnp[n]*vnp[n];

      sum2 += double (n + 1) * hnm[n]*hnm[n] * (pow (r, double(2*n)) - 1.) + double (n - 1) * hnp[n]*hnp[n] * (1. - pow (r, double(-2*n)))
   	    + double (n + 1) * vnm[n]*vnm[n] * (pow (r, double(2*n)) - 1.) + double (n - 1) * vnp[n]*vnp[n] * (1. - pow (r, double(-2*n)));
    }

  delete[] hnp; delete[] hnm; delete[] vnp; delete[] vnm; 
  
  return f1a * log (r) + epsa*epsa * f3a * log (r) - 0.5 * epsa*epsa * f1a * (- sum1 * log (r) + sum2);
}

// ######################
// Function to return f_R
// ######################
double Equilibrium::Getf_R (double r, double w)
{
  if (r >= rc)
    return Getf_R (rc - eps, w) * r*r /rc/rc;
  else
    {
      double L = GetL (r, 0);

      double* hn = new double[Ns+1];
      for (int n = 1; n <= Ns; n++)
	if (r >= 1.)
	  hn[n] = GetHHvac (n, r);
	else
	  hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
      
      double* vn = new double[Ns+1];
      for (int n = 2; n <= Ns; n++)
	if (r >= 1.)
	  vn[n] = GetVVvac (n, r);
	else
	  vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
      
      double* hnp = new double[Ns+1];
      for (int n = 1; n <= Ns; n++)
	if (r >= 1.)
	  hnp[n] = GetHPvac (n, r);
	else
	  hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
      
      double* vnp = new double[Ns+1];
      for (int n = 2; n <= Ns; n++)
	if (r >= 1.)
	  vnp[n] = GetVPvac (n, r);
	else
	  vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
      
      double R = epsa*epsa*epsa * L * cos (w) + epsa*epsa * hn[1];
      
      for (int n = 2; n <= Ns; n++)
	R += epsa*epsa * (hn[n] * cos (double (n - 1) * w) + vn[n] * sin (double (n - 1) * w));
      
      delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
      
      return R;
    }
}

// ######################
// Function to return f_Z
// ######################
double Equilibrium::Getf_Z (double r, double w)
{
  if (r >= rc)
    return Getf_Z (rc - eps, w) * r*r /rc/rc;
  else
    {
      double L = GetL (r, 0);
      
      double* hn = new double[Ns+1];
      for (int n = 1; n <= Ns; n++)
	if (r >= 1.)
	  hn[n] = GetHHvac (n, r);
	else
	  hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
      
      double* vn = new double[Ns+1];
      for (int n = 2; n <= Ns; n++)
	if (r >= 1.)
	  vn[n] = GetVVvac (n, r);
	else
	  vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
      
      double* hnp = new double[Ns+1];
      for (int n = 1; n <= Ns; n++)
	if (r >= 1.)
	  hnp[n] = GetHPvac (n, r);
	else
	  hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
      
      double* vnp = new double[Ns+1];
      for (int n = 2; n <= Ns; n++)
	if (r >= 1.)
	  vnp[n] = GetVPvac (n, r);
	else
	  vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
      
      double Z = - epsa*epsa*epsa * L * sin (w);
      
      for (int n = 2; n <= Ns; n++)
	Z += epsa*epsa * (hn[n] * sin (double (n - 1) * w) - vn[n] * cos (double (n - 1) * w));
      
      delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
      
      return Z;
    }
}

// ##############################
// Function to return |nabla r|^2
// ##############################
double Equilibrium::Getgrr2 (double r, double t)
{
  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hn[n] = GetHHvac (n, r);
    else
      hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vn[n] = GetVVvac (n, r);
    else
      vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    if (r >= 1.)
      hnp[n] = GetHPvac (n, r);
    else
      hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    if (r >= 1.)
      vnp[n] = GetVPvac (n, r);
    else
      vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  double grr2 = 1. + 2. * epsa * hnp[1] * cos (t) + epsa*epsa * (0.75*r*r - hn[1] + hnp[1]*hnp[1] /2.);

  for (int n = 2; n <= Ns; n++)
    {
      grr2 += 2.*epsa * (hnp[n] * cos (double (n) * t) + vnp[n] * sin (double (n) * t))
	+ epsa*epsa * (  (hnp[n]*hnp[n] + double (n*n - 1) * hn[n]*hn[n] /r/r) /2.
		       + (vnp[n]*vnp[n] + double (n*n - 1) * vn[n]*vn[n] /r/r) /2.);
    }

  delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
  
  return grr2;
}

// ########################################
// Function to return theta(omega) function
// ########################################
double Equilibrium::Gettheta (double r, double w, int order)
{
  // .....................
  // Get shaping functions
  // .....................
  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  // ...................
  // Lowest-order result
  // ...................
  double tfun = epsa * r * sin(w) - epsa * hnp[1] * sin(w);

  for (int n = 2; n <= Ns; n++)
    {
      tfun -= epsa * (hnp[n] - double (n - 1) * hn[n] /r) * sin (double (n) * w) /double (n);
      tfun += epsa * (vnp[n] - double (n - 1) * vn[n] /r) * cos (double (n) * w) /double (n);
    }

  // .......................
  // Higher-order correction
  // .......................
  if (order)
    {
      double* sums = new double[Ns+1];
      
      for (int n = 1; n <= Ns; n++)
	{
	  sums[n] = 0.;
	  
	  for (int np = 2; np <= Ns; np++)
	    {
	      if (np + n <= Ns)
		sums[n] +=
		  + (double (np + n - 1) /double (n)) * (hnp[np]   * hn[np+n] + vnp[np]   * vn[np+n]) /r
		  + (double (np     - 1) /double (n)) * (hnp[np+n] * hn[np]   + vnp[np+n] * vn[np])   /r;
	    }
	}
      
      if (Ns >= 2)
	tfun += epsa*epsa * (vn[2] /2. + r*vnp[2] /2. + hnp[1]*vn[2] /r) * cos (w);
      if (Ns >= 3)
	tfun += epsa*epsa * (            r*vnp[3] /4. + hnp[1]*vn[3] /r) * cos (2. * w);
      
      for (int n = 3; n <= Ns; n++)
	{
	  if (n + 1 <= Ns)
	    tfun +=
	      + epsa*epsa * (( - double (n - 2) * (vn[n-1] + vnp[n+1])
			       + r * (vnp[n-1] + vnp[n+1])) /(2. * double (n)) + hnp[1] * vn[n+1] /r) * cos (double (n) * w);
	  else
	    tfun +=
	      + epsa*epsa * ((- double (n - 2) * vn[n-1] + r * vnp[n-1]) /(2. * double (n))) * cos (double (n) * w);
	}
      
      if (Ns >= 2)
	tfun -= epsa*epsa * (hn[2] /2. + r * hnp[2] /2. + hnp[1] * hn[2] /r + sums[1]) * sin (w);
      
      tfun += epsa*epsa * (r*r /4. - r*hnp[1] /4.) * sin (2. * w);
      
      if (Ns >= 3)
	tfun -= epsa*epsa * (r*hnp[3] /4. + hnp[1]*hn[3] /r + sums[2]) * sin (2. * w);
      
      for (int n = 3; n <= Ns; n++)
	{
	  tfun   -= epsa*epsa * ( - (double (n - 2) /double (2*n)) * hn[n-1] + r * hnp[n-1] /double (2*n) + sums[n]) * sin (double (n) * w);
	  
	  if (n + 1 >= Ns)
	    tfun -= epsa*epsa * ( - (double (n - 2) /double (2*n)) * hn[n+1] + r * hnp[n+1] /double (2*n) + hnp[1] * hnp[n+1] /r) * sin (double (n) * w);
	}
      
      delete[] sums; delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
    }
  
  return tfun;
}

// #########################################
// Function to return dtheta/domega function
// #########################################
double Equilibrium::Getdthetadomega (double r, double w, int order)
{
  // .....................
  // Get shaping functions
  // .....................
  double* hn = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    hn[n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
  
  double* vn = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    vn[n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
  
  double* hnp = new double[Ns+1];
  for (int n = 1; n <= Ns; n++)
    hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
  
  double* vnp = new double[Ns+1];
  for (int n = 2; n <= Ns; n++)
    vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
  
  // ...................
  // Lowest-order result
  // ...................
  double tfun = epsa * r * cos(w) - epsa * hnp[1] * cos(w);

  for (int n = 2; n <= Ns; n++)
    {
      tfun -= epsa * (hnp[n] - double (n - 1) * hn[n] /r) * cos (double (n) * w);
      tfun -= epsa * (vnp[n] - double (n - 1) * vn[n] /r) * sin (double (n) * w);
    }

  // .......................
  // Higher-order correction
  // .......................
  if (order)
    {
      double* sums = new double[Ns+1];
      
      for (int n = 1; n <= Ns; n++)
	{
	  sums[n] = 0.;
	  
	  for (int np = 2; np <= Ns; np++)
	    {
	      if (np + n <= Ns)
		sums[n] +=
		  + (double (np + n - 1) /double (n)) * (hnp[np]   * hn[np+n] + vnp[np]   * vn[np+n]) /r
		  + (double (np     - 1) /double (n)) * (hnp[np+n] * hn[np]   + vnp[np+n] * vn[np])   /r;
	    }
	}
      
      if (Ns >= 2)
	tfun -= epsa*epsa * (vn[2] /2. + r*vnp[2] /2. + hnp[1]*vn[2] /r) * sin (w);
      if (Ns >= 3)
	tfun -= epsa*epsa * (            r*vnp[3] /4. + hnp[1]*vn[3] /r) * 2. * sin (2. * w);
      
      for (int n = 3; n <= Ns; n++)
	{
	  if (n + 1 <= Ns)
	    tfun -=
	      + epsa*epsa * (( - double (n - 2) * (vn[n-1] + vnp[n+1])
			       + r * (vnp[n-1] + vnp[n+1])) /(2. * double (n)) + hnp[1] * vn[n+1] /r) * double (n) * sin (double (n) * w);
	  else
	    tfun -=
	      + epsa*epsa * ((- double (n - 2) * vn[n-1] + r * vnp[n-1]) /(2. * double (n))) * double (n) * sin (double (n) * w);
	}
      
      if (Ns >= 2)
	tfun  -= epsa*epsa * (hn[2] /2. + r * hnp[2] /2. + hnp[1] * hn[2] /r + sums[1]) * cos (w);
      
      tfun += epsa*epsa * (r*r /4. - r*hnp[1] /4.) * 2. * cos (2. * w);
      
      if (Ns >= 3)
	tfun -= epsa*epsa * (r*hnp[3] /4. + hnp[1]*hn[3] /r + sums[2]) * 2. * cos (2. * w);
      
      for (int n = 3; n <= Ns; n++)
	{
	  tfun   -= epsa*epsa * ( - (double (n - 2) /double (2*n)) * hn[n-1] + r * hnp[n-1] /double (2*n) + sums[n]) * double (n) * cos (double (n) * w);
	  
	  if (n + 1 >= Ns)
	    tfun -= epsa*epsa * ( - (double (n - 2) /double (2*n)) * hn[n+1] + r * hnp[n+1] /double (2*n) + hnp[1] * hnp[n+1] /r) * double (n) * cos (double (n) * w);
	}
      
      delete[] sums; delete[] hn; delete[] vn; delete[] hnp; delete[] vnp;
    }
  
  return tfun;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Equilibrium::CashKarp45Rhs (double r, double* y, double* dydr)
{
  if (rhs_chooser == 0)
    {
      // ...............................................
      // Right-hand sides for shaping function equations
      // ...............................................

      double* Hn   = new double[Ns+1];
      double* Hnp  = new double[Ns+1];
      double* Hnpp = new double[Ns+1];
      double* Vn   = new double[Ns+1];
      double* Vnp  = new double[Ns+1];
      double* Vnpp = new double[Ns+1];
      
      double g2 = y[0];

      Hn [1] = y[1];
      Hnp[1] = y[2];
      int j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  Hn [n] = y[j]; j++;
	  Hnp[n] = y[j]; j++;
	  Vn [n] = y[j]; j++;
	  Vnp[n] = y[j]; j++;
	}

      double f1  = Getf1 (r);
      double f1p = Getf1p(r);
      double p2p = Getp2p(r);
      
      double facf = 2.*f1p/f1 - 1./r;
      double facp = 2.*r*r*r*p2p/f1/f1;
      
      double g2p = - p2p - f1*f1p/r/r;

      Hnpp[1] = - facf * Hnp[1] - 1. + facp;

      for (int n = 2; n <= Ns; n++)
	{
	  Hnpp[n] = - facf * Hnp[n] + double (n*n - 1) * Hn[n]/r/r;
	  Vnpp[n] = - facf * Vnp[n] + double (n*n - 1) * Vn[n]/r/r;
  	}
      
      dydr[0] = g2p;
      dydr[1] = Hnp [1];
      dydr[2] = Hnpp[1];
      j = 3;
      for (int n = 2; n <= Ns; n++)
	{
	  dydr[j] = Hnp [n]; j++;
	  dydr[j] = Hnpp[n]; j++;
	  dydr[j] = Vnp [n]; j++;
	  dydr[j] = Vnpp[n]; j++;
	}

      delete[] Hn; delete[] Hnp; delete[] Hnpp; delete[] Vn; delete[] Vnp; delete[] Vnpp;
     }
  else if (rhs_chooser == 1)
    {
      // ................................
      // Right-hand side for f13 equation
      // ................................

      double* Hn  = new double[Ns+1];
      double* Hnp = new double[Ns+1];
      double* Vn  = new double[Ns+1];
      double* Vnp = new double[Ns+1];

      double f1  = Getf1 (r);
      double f1p = Getf1p(r);
      double p2p = Getp2p(r);
      
      double g2 = gsl_spline_eval (g2spline, r, g2acc);

      Hn [1] = gsl_spline_eval (HHspline[1], r, HHacc[1]);
      Hnp[1] = gsl_spline_eval (HPspline[1], r, HPacc[1]);
      for (int n = 2; n <= Ns; n++)
	{
	  Hn [n] = gsl_spline_eval (HHspline[n], r, HHacc[n]);
	  Hnp[n] = gsl_spline_eval (HPspline[n], r, HPacc[n]);
	  Vn [n] = gsl_spline_eval (VVspline[n], r, VVacc[n]);
	  Vnp[n] = gsl_spline_eval (VPspline[n], r, VPacc[n]);
	}

      double fac1 = 0., fac2 = 0.;
      for (int n = 2; n <= Ns; n++)
	{
	  fac1 += Hnp[n]*Hnp[n] + 2. * double (n*n - 1) * Hnp[n]*Hn[n]/r - double (n*n - 1) * Hn[n]*Hn[n]/r/r;
	  fac1 += Vnp[n]*Vnp[n] + 2. * double (n*n - 1) * Vnp[n]*Vn[n]/r - double (n*n - 1) * Vn[n]*Vn[n]/r/r;

	  fac2 += (3.*Hnp[n]*Hnp[n] - double (n*n - 1) * Hn[n]*Hn[n]/r/r) /2.;
	  fac2 += (3.*Vnp[n]*Vnp[n] - double (n*n - 1) * Vn[n]*Vn[n]/r/r) /2.;
	}

      double f3p =
	- y[0]*f1p/f1
	- f1 * (3.*r*r/2. - 2.*r*Hnp[1] + Hnp[1]*Hnp[1] + fac1)/r
	+ f1p * (g2 - 3.*r*r/4. + Hn[1] + 3.*Hnp[1]*Hnp[1]/2. + fac2)
	+ r*r*p2p * (g2 + r*r/2. - 3.*r*Hnp[1] - 2.*Hn[1]) /f1;
 
      dydr[0] = f3p;
      dydr[1] = (f1 + epsa*epsa * y[0]) /r;
 
      delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;
    }
  else if (rhs_chooser == 2)
    {
      // ..........................................
      // Right-hand sides for li and beta equations
      // ..........................................

      double p2  = Getp2 (r);
      double g2  = gsl_spline_eval (g2spline,  r, g2acc);
      double f   = gsl_spline_eval (fspline,   r, facc);
      double gr2 = gsl_spline_eval (gr2spline, r, gr2acc);
      double R2  = gsl_spline_eval (R2spline,  r, R2acc);

      dydr[0] = f*f * gr2 /r;
      dydr[1] = r * p2 * R2;
      dydr[2] = r * (1. + eps*eps * g2) * (1. + eps*eps * g2);
      dydr[3] = r * R2;
    }
  else if (rhs_chooser == 3)
    {
      // ........................................
      // Right-hand side for boundary calculation
      // ........................................

      double rb    = 1.;
      double R     = GetR    (rb, r, 1);
      double dRdr  = GetdRdr (rb, r, 1);
      double dRdw  = GetdRdw (rb, r, 1);
      double dZdr  = GetdZdr (rb, r, 1);
      double dZdw  = GetdZdw (rb, r, 1);

      dydr[0] = (dRdw * dZdr - dRdr * dZdw) /R;
    }
  else if (rhs_chooser == 4)
    {
      // ....................................
      // Right-hand side for wall calculation
      // ....................................

      double rw    = bw;
      double R     = GetR    (rw, r, 1);
      double dRdr  = GetdRdr (rw, r, 1);
      double dRdw  = GetdRdw (rw, r, 1);
      double dZdr  = GetdZdr (rw, r, 1);
      double dZdw  = GetdZdw (rw, r, 1);

      dydr[0] = (dRdw * dZdr - dRdr * dZdw) /R;
    }
 }
