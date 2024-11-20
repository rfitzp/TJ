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
  string JSONFilename = "../Inputs/Equilibrium.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

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

  for (const auto& number : JSONData["Hna"])
    {
      Hna.push_back (number.get<double> ());
    }
  for (const auto& number : JSONData["Vna"])
    {
      Vna.push_back (number.get<double> ());
    }

  JSONFilename = "../Inputs/Layer.json";
  JSONData     = ReadJSONFile (JSONFilename);

  B0 = JSONData["B0"].get<double> ();
  R0 = JSONData["R0"].get<double> ();

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

// #################################################
// Function to override qc value read from JSON file
// #################################################
void Equilibrium::Setqc (double _qc)
{
  qc = _qc;
}
// #################################################
// Function to override qa value read from JSON file
// #################################################
void Equilibrium::Setqa (double _qa)
{
  qa = _qa;
}

// ###################################################
// Function to override epsa value read from JSON file
// ###################################################
void Equilibrium::Setepsa (double _epsa)
{
  epsa = _epsa;
}

// #################################################
// Function to override pc value read from JSON file
// #################################################
void Equilibrium::Setpc (double _pc)
{
  pc = _pc;
}

// #####################################################
// Function to override Hna[0] value read from JSON file
// #####################################################
void Equilibrium::SetH2 (double _H2)
{
  Hna[0] = _H2;
}

// #####################################################
// Function to override Vna[0] value read from JSON file
// #####################################################
void Equilibrium::SetV2 (double _V2)
{
  Vna[0] = _V2;
}
// #####################################################
// Function to override Hna[1] value read from JSON file
// #####################################################
void Equilibrium::SetH3 (double _H3)
{
  Hna[1] = _H3;
}

// #####################################################
// Function to override Vna[1] value read from JSON file
// #####################################################
void Equilibrium::SetV3 (double _V3)
{
  Vna[1] = _V3;
}
// #####################################################
// Function to override Hna[2] value read from JSON file
// #####################################################
void Equilibrium::SetH4 (double _H4)
{
  Hna[2] = _H4;
}

// #####################################################
// Function to override Vna[2] value read from JSON file
// #####################################################
void Equilibrium::SetV4 (double _V4)
{
  Vna[2] = _V4;
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
  printf ("qc = %10.3e qa = %10.3e epsa = %10.3e pc = %10.3e mu = %10.3e Ns = %3d Nr = %3d Nf = %3d Nw = %3d\n",
	  qc, qa, epsa, pc, mu, Ns, Nr, Nf, Nw);

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
  Lfunc = new double[Nr+1];

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
  twall    = new double[Nw+1];
  wwall0   = new double[Nw+1];
  twall0   = new double[Nw+1];
  wwall    = new double[Nw+1];
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
      Rhs (rr[i], y1, dy1dr);
      
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
	  gr2 += (HPfunc(n, i) * HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i) * HHfunc(n, i) /r/r)/2.;
	  gr2 += (VPfunc(n, i) * VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i) * VVfunc(n, i) /r/r)/2.;
	  ir2 += (HPfunc(n, i) * HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i) * HHfunc(n, i) /r/r)/2.;
	  ir2 += (VPfunc(n, i) * VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i) * VVfunc(n, i) /r/r)/2.;
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
      Jt[i] = gsl_spline_eval_deriv (Itspline, rr[i], Itacc);
      Jp[i] = gsl_spline_eval_deriv (Ipspline, rr[i], Ipacc);
      DI[i] = GetDI (rr[i]);
      DR[i] = GetDR (rr[i]);
    }
  DI[0] = DI[1];
  DR[0] = DR[1];

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

  printf ("qc = %10.3e q0a   = %10.3e q2a   = %10.3e Ip    = %10.3e It = %10.3e\n",
	  q2[0], q0[Nr], q2[Nr], Ip[Nr], It[Nr]);
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
 
      for (int n = 2; n <= Ns; n++)
	{
	  sum1 += (                    3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		       - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]) /2.;
	  sum2 +=                            HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i)
	          + 2. * double (n*n - 1) * (HPfunc(n, i) * HHfunc(n, i) + VPfunc(n, i) * VVfunc(n, i)) /rr[i]
	          -      double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i];
	  sum3 += (                    3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		       - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]) /2.;
	}

      S1[i] = sum1;
      S2[i] = sum2;
      S3[i] = sum3;
    }
  S1[0] = 0.;
  S2[0] = S2[1];
  S3[0] = 0.;

  // ...........................
  // Calculate profile functions
  // ...........................
  for (int i = 0; i <= Nr; i++)
    {
      P1[i]  = (2. - s[i]) /q2[i];
      P1a[i] = (2. - s0[i]) /q0[i];
      P2[i]  = (- 3. * s[i] + 2. * s[i]*s[i] - s2[i]) /q2[i];
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
	  
	  double w, wold = t;
	  for (int i = 0; i < 10; i++)
	    {
	      w    = t - Gettheta (rf, wold, 1);
	      wold = w;
	    }
	  
	  double R = GetR (rf, w, 1);
	  double Z = GetZ (rf, w, 1);
	  
	  RR    (i-1, j) = R;
	  ZZ    (i-1, j) = Z;
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
	  
	  RRw (i-1, j) = R;
	  ZZw (i-1, j) = Z;
	}
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

  delete[] y3; delete[] err3;
  
  double rbb = tbound0[Nw] /(2.*M_PI);
  printf ("rb = %10.3e\n", rbb);
  
  for (int j = 0; j < Nw; j++)
    tbound0[j] /= rbb;

  tbound0[Nw] = 2.*M_PI;
  
  // Interpolate preliminary boundary data
  gsl_spline_init (wspline, tbound0, wbound0, Nw+1);
  
  // Calculate final boundary data
  for (int j = 0; j <= Nw; j++)
    tbound[j] = double (j) * 2.*M_PI /double (Nw);
  
  double rb = 1.;
  for (int j = 0; j <= Nw; j++)
    {
      double t = tbound[j];
      double w = gsl_spline_eval (wspline, t, wacc);
      
      wbound[j] = w;
      Rbound[j] = GetR    (rb, w, 1);
      Zbound[j] = GetZ    (rb, w, 1);
      R2b   [j] = GetR2   (rb, t);
      grr2b [j] = Getgrr2 (rb, t);
    }
  
  gsl_spline_init (Rspline, tbound, Rbound, Nw+1);
  gsl_spline_init (Zspline, tbound, Zbound, Nw+1);
  
  for (int j = 0; j <= Nw; j++)
    {
      dRdtheta[j] = gsl_spline_eval_deriv (Rspline, tbound[j], Racc);
      dZdtheta[j] = gsl_spline_eval_deriv (Zspline, tbound[j], Zacc);
    }
   
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
  delete[] s0;  
 
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

  delete[] Rwall; delete[] Zwall; delete[] twall; delete[] wwall0;   delete[] twall0;
  delete[] wwall; delete[] R2w;   delete[] grr2w; delete[] dRdthetw; delete[] dZdthetw;

  delete[] Psi; delete[] PsiN; delete[] Tf;    delete[] mu0P;
  delete[] DI;   delete[] DR;  delete[] Lfunc;  
 } 


// ########################
// Function to return f1(r)
// ########################
double Equilibrium::Getf1 (double r)
{
  if (r < 0.1)
    return r*r * (1. - (nu-1.)*r*r/2. + (nu-1.)*(nu-2.)*r*r*r*r/6. - (nu-1.)*(nu-2.)*(nu-3.)*r*r*r*r*r*r/24.)/qc;
  else
    return (1. - pow (1. - r*r, nu)) /nu/qc;
}

// #########################
// Function to return f1'(r)
// #########################
double Equilibrium::Getf1p (double r)
{
  if (r < 0.1)
    return 2.*r * (1. - (nu-1.)*r*r + (nu-1.)*(nu-2.)*r*r*r*r/2. - (nu-1.)*(nu-2.)*(nu-3.)*r*r*r*r*r*r/6.)/qc;
  else
    return 2. * r * pow (1. - r*r, nu-1.) /qc;
}

// ########################
// Function to return p2(r)
// ########################
double Equilibrium::Getp2 (double r)
{
  return pc * pow (1. - r*r, mu);
}

// #########################
// Function to return p2'(r)
// #########################
double Equilibrium::Getp2p (double r)
{
  return - 2. * pc * mu * r * pow (1. - r*r, mu-1.);
}

// ##########################
// Function to return p2''(r)
// ##########################
double Equilibrium::Getp2pp (double r)
{
  return - 2. * pc * mu * pow (1. - r*r, mu-1.)
    + 4. * pc * mu * (mu-1.) * r*r * pow (1. - r*r, mu-2.);
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
	tfun  -= epsa*epsa * (hn[2] /2. + r * hnp[2] /2. + hnp[1] * hn[2] /r + sums[1]) * sin (w);
      
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

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Equilibrium::Rhs (double r, double* y, double* dydr)
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
 }

// ###################################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using Cash-Karp adaptive step-length fourth-order/fifth-order Runge-Kutta scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... safety factor
//     T       ... step-length cannot change by more than this factor from step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step-length
//     h_max   ... maximum allowable step-length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step while attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
// ####################################################################################
void Equilibrium::CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
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

  // Take Cash-Karp RK4/RK5 step 
  CashKarp45Fixed (neqns, x, y, Err, h);

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

  // Calculate new step-length
  double h_est;
  if (acc >= t_err)
    h_est = S * h * pow (fabs (acc / t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc / t_err), 0.25);

  // Prevent step-length from changing by more than factor T
  if (h_est / h > T)
    h *= T;
  else if (h_est / h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Prevent step-length from falling below h_min
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
      CashKarp45Adaptive (neqns, x, y, h, t_err, acc, S, T, rept, 
			  maxrept, h_min, h_max, flag, diag, file);
    }

  delete[] y0; 
  delete[] Err;
}

// #####################################################################3
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length Cash-Karp fourth-order/fifth-order Runge-Kutta
// scheme
//
//     neqns   ... number of equations
//     x       ... independent variable
//     y       ... array of dependent variables 
//     err     ... array of errors
//     h       ... step-length
//     
// ######################################################################
void Equilibrium::CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h)
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

// ########################################
// Function to strip comments from a string
// ########################################
string Equilibrium::stripComments (const string& input)
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
json Equilibrium::ReadJSONFile (const string& filename)
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
FILE* Equilibrium::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Equilibrium::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* Equilibrium::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("Equilibrium::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ################################################################
// Function to check that directory exists, and create it otherwise
// ################################################################
bool Equilibrium::CreateDirectory (const char* path)
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
