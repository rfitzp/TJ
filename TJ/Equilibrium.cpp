// Equilibrium.cpp

#include "LightEquilibrium.h"
#include "Equilibrium.h"

// ###########
// Constructor
// ###########
Equilibrium::Equilibrium ()
{
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
  string JSONFilename = "Inputs/Equilibrium.json";
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
  if (epsa < 0.)
    {
      printf ("Equilibrium:: Error - epsa cannot be less than unity\n");
      exit (1);
    }
  if (eps < 0.)
    {
      printf ("Equilibrium:: Error - eps cannot be less than unity\n");
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
  LightEquilibrium lightequilibrium;
      
  nu = lightequilibrium.GetNu (qa);
}

// #################################################
// Function to override pc value read from JSON file
// #################################################
void Equilibrium::Setpc (double _pc)
{
  pc = _pc;
}

// #####################################
// Function to solve equilibrium problem
// #####################################
void Equilibrium::Solve ()
{
  // ...............
  // Allocate memory
  // ...............
  rr   = new double[Nr+1];
  p2   = new double[Nr+1];
  f1   = new double[Nr+1];
  f3   = new double[Nr+1];
  g2   = new double[Nr+1];
  q0   = new double[Nr+1];
  q2   = new double[Nr+1];
  It   = new double[Nr+1];
  Ip   = new double[Nr+1];
  Jt   = new double[Nr+1];
  Jp   = new double[Nr+1];
  pp   = new double[Nr+1];
  ppp  = new double[Nr+1];
  qq   = new double[Nr+1];
  qqq  = new double[Nr+1];
  s    = new double[Nr+1];
  s2   = new double[Nr+1];
  S1   = new double[Nr+1];
  S2   = new double[Nr+1];
  P1   = new double[Nr+1];
  P2   = new double[Nr+1];
  P3   = new double[Nr+1];
  P3a  = new double[Nr+1];
  ff   = new double[Nr+1];
  ggr2 = new double[Nr+1];
  RR2  = new double[Nr+1];
  Psi  = new double[Nr+1];

  HHfunc.resize (Ns+1, Nr+1);
  VVfunc.resize (Ns+1, Nr+1);
  HPfunc.resize (Ns+1, Nr+1);
  VPfunc.resize (Ns+1, Nr+1);
  Lfunc = new double[Nr+1];

  Itspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Ipspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  g2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  fspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  q2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  gr2spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  R2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
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
  
  printf ("\n");
  printf ("Class EQUILIBRIUM::\n");
  printf ("Calculation parameters:\n");
  printf ("qc  = %10.3e nu = %10.3e pc = %10.3e mu = %10.3e epsa = %10.3e Ns = %3d Nr = %3d Nf = %3d Nw = %3d\n",
	  qc, nu, pc, mu, epsa, Ns, Nr, Nf, Nw);
  
  // .........................
  // Rescale shaping functions
  // .........................
  int    nn   = 1;
  double zero = 0.;
  printf ("n = %3d:  Hna = %10.3e  Vna = %10.3e\n", nn, HHfunc(1, Nr), zero);

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
      q2[i]  = rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) /(f1[i] + epsa*epsa*f3[i]);

      double ff1 = f1[i];
      double ff3 = f3[i];
      double f1p = Getf1p (r);
      double p2p = pp[i];
      double gg2 = g2[i];
      double g2p = - p2p - ff1*f1p/r/r;
      double f3p = dy1dr[0];

      qq[i] = 2.*rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) /(f1[i] + epsa*epsa*f3[i])
	+ rr[i]*rr[i]*rr[i] * (epsa*epsa*g2p) /(f1[i] + epsa*epsa*f3[i])
	- rr[i]*rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) * (f1p + epsa*epsa*f3p)
	/(f1[i] + epsa*epsa*f3[i]) /(f1[i] + epsa*epsa*f3[i]);
  
      double gr2 = 3.*rr[i]*rr[i]/4. - HHfunc(1, i) + HPfunc(1, i) * HPfunc(1, i) /2.;
      for (int n = 2; n <= Ns; n++)
	{
	  gr2 += (HPfunc(n, i) * HPfunc(n, i) + double (n*n - 1) * HHfunc(n, i) * HHfunc(n, i) /r/r)/2.;
	  gr2 += (VPfunc(n, i) * VPfunc(n, i) + double (n*n - 1) * VVfunc(n, i) * VVfunc(n, i) /r/r)/2.;
	}

      double R2 = rr[i]*rr[i]/2. - rr[i]*HPfunc(1, i) - 2.*HHfunc(1, i);

      ff  [i] = ff1 + epsa*epsa * ff3;
      ggr2[i] = 1.  + epsa*epsa * gr2;
      RR2 [i] = 1.  - epsa*epsa * R2;
      It  [i] =   2.*M_PI * (f1[i] * (1. + epsa*epsa * gr2) + epsa*epsa * f3[i]);
      Ip  [i] = - 2.*M_PI * g2[i];
    }
  q2[0] = qc * (1. + epsa*epsa * (H2c*H2c + V2c*V2c));

  delete[] y1; delete[] dy1dr; delete[] err1;

  // ........................
  // Calculate magnetic shear
  // ........................
  for (int i = 0; i <= Nr; i++)
    {
      s[i] = qq[i] /q2[i];
    }
  
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
   double sa = 2. + epsa*epsa * sum * q0[Nr] /q2[Nr];
   
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
  gsl_spline_init (Itspline,  rr, It,   Nr+1);
  gsl_spline_init (Ipspline,  rr, Ip,   Nr+1);

  for (int i = 0; i <= Nr; i++)
    {
      Jt[i] = gsl_spline_eval_deriv (Itspline, rr[i], Itacc);
      Jp[i] = gsl_spline_eval_deriv (Ipspline, rr[i], Ipacc);
    }

  // ............................
  // Calculate li and beta values
  // ............................
  double* y2   = new double[3];
  double* err2 = new double[3];
  rhs_chooser  = 2;

  r     = eps;
  h     = h0;
  count = 0;
  y2[0] = 0.;
  y2[1] = 0.;
  y2[2] = 0.;

  do
    {
      CashKarp45Adaptive (3, r, y2, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < 1. - h);
  CashKarp45Fixed (3, r, y2, err2, 1. - r);

  li    = 2.*y2[0] /ff[Nr]/ff[Nr] /ggr2[Nr] /ggr2[Nr];
  betat = 2. * epsa*epsa * y2[1] /y2[2];
  betap = 2. * y2[1] /y2[0];
  betaN = 20. * betat /epsa /ff[Nr] /ggr2[Nr];

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
      double sum2 = - (1. - s[i]) * (HPfunc(1, i) * HHfunc(1, i) /rr[i] + (1./3.) * HHfunc(1, i) * HHfunc(1, i) /rr[i]/rr[i]);
 
      for (int n = 2; n <= Ns; n++)
	{
	  sum1 += (                  3. * (HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i))
		     - double (n*n - 1) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i])/2.;
	  sum2 += 2. * double (n*n - 1) * (+ HPfunc(n, i) * HPfunc(n, i) + VPfunc(n, i) * VPfunc(n, i)
					- (11./3.) * (HPfunc(n, i) * HHfunc(n, i) + VPfunc(n, i) * VVfunc(n, i)) /rr[i]
					+ double (n*n) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i)*VVfunc(n, i)) /rr[i]/rr[i])
	    - (1. - s[i]) * (            (HPfunc(n, i) * HHfunc(n, i) + VPfunc(n, i) * VVfunc(n, i)) /rr[i]
			     + (1./3.) * (HHfunc(n, i) * HHfunc(n, i) + VVfunc(n, i) * VVfunc(n, i)) /rr[i]/rr[i]);
	}

      S1[i] = sum1;
      S2[i] = sum2;
    }
  S1[0] = 0.;
  S2[0] = S2[1];

  // ...........................
  // Calculate profile functions
  // ...........................
  for (int i = 0; i <= Nr; i++)
    {
      P1[i]  = (2. - s[i]) /q2[i];
      P2[i]  = (- 3. * s[i] + 2. * s[i]*s[i] - s2[i]) /q2[i];
      P3a[i] =
        + rr[i]*rr[i]   * (2. - s[i]) /q2[i]/q2[i]/q2[i]
	+ (s[i] /q2[i]) * (3. * rr[i]*rr[i] /4. - HHfunc(1, i) - S1[i])
	- (2.   /q2[i]) * (3. * rr[i]*rr[i] /2. - HHfunc(1, i) - rr[i] * HPfunc(1, i) - 2. * S1[i] /3.);
    }
  P1[0] = P1[1];

  for (int i = 1; i < Nr; i++)
    {
      P3[i] = 2. * rr[i] * pp[i] * (2. - s[i]) + q2[i] * rr[i] * (P3a[i+1] - P3a[i-1]) /2./hh - S2[i];
    }
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
  int serr = system ("./write_efit_");
  if (serr != 0)
    printf ("TJ:: Warning: Error runing WriteEFIT\n");
  
  // ........
  // Clean up
  // ........
  printf ("Cleaning up:\n");
  
  delete[] rr;  delete[] p2; delete[] f1;  delete[] f3; delete[] g2;
  delete[] q0;  delete[] q2; delete[] It;  delete[] Ip; delete[] Jt;
  delete[] Jp;  delete[] pp; delete[] ppp; delete[] qq; delete[] qqq;
  delete[] s;   delete[] s2; delete[] S1;  delete[] S2; delete[] P1;
  delete[] P2;  delete[] P3; delete[] P3a; delete[] ff; delete[] ggr2;
  delete[] RR2;
 
  gsl_spline_free (Itspline);
  gsl_spline_free (Ipspline);
  gsl_spline_free (g2spline);
  gsl_spline_free (fspline);
  gsl_spline_free (q2spline);
  gsl_spline_free (gr2spline);
  gsl_spline_free (R2spline);
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

  delete[] Rbound;   delete[] Zbound; delete[] tbound; delete[] wbound0; delete[] tbound0;
  delete[] wbound;   delete[] R2b;    delete[] grr2b;  delete[] Lfunc;   delete[] dRdtheta; 
  delete[] dZdtheta; delete[] Psi;
 }

// #################################################
// Function to write equilibrium data to netcdf file
// #################################################
void Equilibrium::WriteNetcdf (double sa)
{
  printf ("Writing equilibrium data to netcdf file Outputs/Equilibrium/Equilibrium.nc:\n");

  double  para[2], Input[15], Beta[4];

  double* Hna   = new double[Ns+1];
  double* Vna   = new double[Ns+1];
  double* npol  = new double[Ns+1];

  Input[0]  = qc;
  Input[1]  = qa;
  Input[2]  = nu;
  Input[3]  = pc;
  Input[4]  = mu;
  Input[5]  = epsa;
  Input[6]  = eps;
  Input[7]  = double (Ns);
  Input[8]  = double (Nr);
  Input[9]  = double (Nf);
  Input[10] = double (Nw);
  Input[11] = acc;
  Input[12] = h0;
  Input[13] = hmin;
  Input[14] = hmax;
  
  para[0] = epsa;
  para[1] = sa;

  Beta[0] = li;
  Beta[1] = betat;
  Beta[2] = betap;
  Beta[3] = betaN;
    
  Hna[0] = 0.;
  Vna[0] = 0.;
  Hna[1] = HHfunc(1, Nr);
  Vna[1] = 0.;
  for (int n = 2; n <= Ns; n++)
    {
      Hna[n] = HHfunc(n, Nr);
      Vna[n] = VVfunc(n, Nr);
    }

  for (int n = 0; n <= Ns; n++)
    npol[n] = double (n);
 
  try
    {
      NcFile dataFile ("Outputs/Equilibrium/Equilibrium.nc", NcFile::replace);

      NcDim i_d = dataFile.addDim ("Ni", 15);
      NcDim p_d = dataFile.addDim ("Np", 2);
      NcDim b_d = dataFile.addDim ("Nb", 4);
      NcDim r_d = dataFile.addDim ("Nr", Nr+1);
      NcDim s_d = dataFile.addDim ("Ns", Ns+1);
      NcDim f_d = dataFile.addDim ("Nf", Nf);
      NcDim w_d = dataFile.addDim ("Nw", Nw+1);
 
      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> flux_d;
      flux_d.push_back (f_d);
      flux_d.push_back (w_d);

      NcVar i_x = dataFile.addVar ("InputParameters", ncDouble, i_d);
      i_x.putVar (Input);
      NcVar p_x = dataFile.addVar ("para",            ncDouble, p_d);
      p_x.putVar (para);
      NcVar b_x = dataFile.addVar ("beta",            ncDouble, b_d);
      b_x.putVar (Beta);

      NcVar r_x   = dataFile.addVar ("r",    ncDouble, r_d);
      r_x.putVar (rr);
      NcVar g2_x  = dataFile.addVar ("g_2",  ncDouble, r_d);
      g2_x.putVar (g2);
      NcVar p2_x  = dataFile.addVar ("p_2",  ncDouble, r_d);
      p2_x.putVar (p2);
      NcVar pp_x  = dataFile.addVar ("pp",   ncDouble, r_d);
      pp_x.putVar (pp);
      NcVar ppp_x = dataFile.addVar ("ppp",  ncDouble, r_d);
      ppp_x.putVar (ppp);
      NcVar f1_x  = dataFile.addVar ("f_1",  ncDouble, r_d);
      f1_x.putVar (f1);
      NcVar f3_x  = dataFile.addVar ("f_3",  ncDouble, r_d);
      f3_x.putVar (f3);
      NcVar q0_x  = dataFile.addVar ("q_0",  ncDouble, r_d);
      q0_x.putVar (q0);
      NcVar q2_x  = dataFile.addVar ("q_2",  ncDouble, r_d);
      q2_x.putVar (q2);
      NcVar It_x  = dataFile.addVar ("I_t",  ncDouble, r_d);
      It_x.putVar (It);
      NcVar Ip_x  = dataFile.addVar ("I_p",  ncDouble, r_d);
      Ip_x.putVar (Ip);
      NcVar Jt_x  = dataFile.addVar ("J_t",  ncDouble, r_d);
      Jt_x.putVar (Jt);
      NcVar Jp_x  = dataFile.addVar ("J_p",  ncDouble, r_d);
      Jp_x.putVar (Jp);
      NcVar q_x   = dataFile.addVar ("q",    ncDouble, r_d);
      q_x.putVar (q2);
      NcVar qq_x  = dataFile.addVar ("qq",   ncDouble, r_d);
      qq_x.putVar (qq);
      NcVar qqq_x = dataFile.addVar ("qqq",  ncDouble, r_d);
      qqq_x.putVar (qqq);
      NcVar s_x   = dataFile.addVar ("s",    ncDouble, r_d);
      s_x.putVar (s);
      NcVar s2_x  = dataFile.addVar ("s2",   ncDouble, r_d);
      s2_x.putVar (s2);
      NcVar S1_x  = dataFile.addVar ("S1",   ncDouble, r_d);
      S1_x.putVar (S1);
      NcVar S2_x  = dataFile.addVar ("S2",   ncDouble, r_d);
      S2_x.putVar (S2);
      NcVar P1_x  = dataFile.addVar ("P1",   ncDouble, r_d);
      P1_x.putVar (P1);
      NcVar P2_x  = dataFile.addVar ("P2",   ncDouble, r_d);
      P2_x.putVar (P2);
      NcVar P3_x  = dataFile.addVar ("P3",   ncDouble, r_d);
      P3_x.putVar (P3);
      NcVar P3a_x = dataFile.addVar ("P3a",  ncDouble, r_d);
      P3a_x.putVar (P3a);
 
      NcVar Hn_x  = dataFile.addVar ("Hn",  ncDouble, shape_d);
      Hn_x.putVar (HHfunc.data());
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (HPfunc.data());
      NcVar Vn_x  = dataFile.addVar ("Vn",  ncDouble, shape_d);
      Vn_x.putVar (VVfunc.data());
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (VPfunc.data());

      NcVar R_x  = dataFile.addVar ("R",     ncDouble, flux_d);
      R_x.putVar (RR.data());
      NcVar Z_x  = dataFile.addVar ("Z",     ncDouble, flux_d);
      Z_x.putVar (ZZ.data());
      NcVar Rw_x = dataFile.addVar ("Rw",    ncDouble, flux_d);
      Rw_x.putVar (RRw.data());
      NcVar Zw_x = dataFile.addVar ("Zw",    ncDouble, flux_d);
      Zw_x.putVar (ZZw.data());
      NcVar rr_x = dataFile.addVar ("rr",    ncDouble, flux_d);
      rr_x.putVar (rvals.data());
      NcVar t_x  = dataFile.addVar ("theta", ncDouble, flux_d);
      t_x.putVar (thvals.data());
      NcVar w_x  = dataFile.addVar ("omega", ncDouble, flux_d);
      w_x.putVar (wvals.data());

      NcVar n_x   = dataFile.addVar ("n",   ncDouble, s_d);
      n_x.putVar (npol);
      NcVar Hna_x = dataFile.addVar ("Hna", ncDouble, s_d);
      Hna_x.putVar (Hna);
      NcVar Vna_x = dataFile.addVar ("Vna", ncDouble, s_d);
      Vna_x.putVar (Vna);

      NcVar Rbound_x   = dataFile.addVar ("Rbound",    ncDouble, w_d);
      Rbound_x.putVar (Rbound);
      NcVar Zbound_x   = dataFile.addVar ("Zbound",    ncDouble, w_d);
      Zbound_x.putVar (Zbound);
      NcVar tbound_x   = dataFile.addVar ("tbound",    ncDouble, w_d);
      tbound_x.putVar (tbound);
      NcVar wbound_x   = dataFile.addVar ("wbound",    ncDouble, w_d);
      wbound_x.putVar (wbound);
      NcVar R2b_x      = dataFile.addVar ("R2bound",   ncDouble, w_d);
      R2b_x.putVar (R2b);
      NcVar grr2b_x    = dataFile.addVar ("grr2bound", ncDouble, w_d);
      grr2b_x.putVar (grr2b);
      NcVar dRdtheta_x = dataFile.addVar ("dRdtheta",  ncDouble, w_d);
      dRdtheta_x.putVar (dRdtheta);
      NcVar dZdtheta_x = dataFile.addVar ("dZdtheta",  ncDouble, w_d);
      dZdtheta_x.putVar (dZdtheta);
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Equilibrium/Equilbrium.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] npol; delete[] Hna; delete[] Vna;
}

// ###############################
// Function to calculate EFIT data
// ###############################
void Equilibrium::CalculateEFIT ()
{
  // ----------------------
  // Set physical constants
  // ----------------------
  double mu_0 = 4.*M_PI*1.e-7;

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  {
    string JSONFilename = "Inputs/Equilibrium.json";
    json   JSONData     = ReadJSONFile (JSONFilename);

    EFIT  = JSONData["EFIT"].get<int>     ();
    NRBOX = JSONData["NRBOX"].get<int>    ();
    NZBOX = JSONData["NZBOX"].get<int>    ();
    rc    = JSONData["rc"]   .get<double> ();
  }
  if (!EFIT)
    return;

  printf ("Calculating EFIT data:\n");
  printf ("NRBOX = %4d        NZBOX = %4d       rc = %10.3e\n", NRBOX, NZBOX, rc);

  {
    string JSONFilename = "Inputs/TJ.json";
    json   JSONData     = ReadJSONFile (JSONFilename);
    
    B0EXP = JSONData["B0"].get<double> ();
    R0EXP = JSONData["R0"].get<double> ();
  }

  printf ("B0EXP = %10.3e  R0EXP = %10.3e\n", B0EXP, R0EXP);

  RAXIS   = R0EXP;
  ZAXIS   = 0.;
  CURRENT = epsa*epsa * (B0EXP * R0EXP /mu_0) * It[Nr];

  // ---------------
  // Allocate memory
  // ---------------
  PSI  = new double[NRBOX];
  PSIN = new double[NRBOX];
  rPSI = new double[NRBOX];
  T    = new double[NRBOX];
  TTp  = new double[NRBOX];
  P    = new double[NRBOX];
  Pp   = new double[NRBOX];
  Q    = new double[NRBOX];

  RBOUND   = new double[Nw+1];
  ZBOUND   = new double[Nw+1];
  RLIMITER = new double[5];
  ZLIMITER = new double[5];

  RGRID = new double[NRBOX];
  ZGRID = new double[NZBOX];
  PSIRZ = new double[NRBOX*NZBOX];
  rRZ   = new double[NRBOX*NZBOX];
  wRZ   = new double[NRBOX*NZBOX];

  rPsispline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  PSIrspline = gsl_spline_alloc (gsl_interp_cspline, NRBOX);
  rPsiacc    = gsl_interp_accel_alloc ();
  PSIracc    = gsl_interp_accel_alloc ();

  // ------------------------
  // Set up Psi and PSI grids
  // ------------------------
  for (int i = 0; i <= Nr; i++)
    Psi[i] *= epsa*epsa * B0EXP * R0EXP*R0EXP;

  for (int i = 0; i <= Nr; i++)
    Psi[i] -= Psi[Nr];

  PSIAXIS  = Psi[0];
  PSIBOUND = Psi[Nr];

  for (int i = 0; i < NRBOX; i++)
    PSI[i] = PSIAXIS * (double (NRBOX - 1 - i) /double (NRBOX - 1));

  for (int i = 0; i < NRBOX; i++)
    PSIN[i] = (PSI[i] - PSIAXIS) /(PSIBOUND - PSIAXIS);

  // ---------------------------
  // Interpolate r onto Psi grid
  // ---------------------------
  gsl_spline_init (rPsispline, Psi, rr, Nr+1);

  // ----------------
  // Set up rPSI grid
  // ----------------
  rPSI[0]       = 0.;
  rPSI[NRBOX-1] = 1.;
  for (int i = 1; i < NRBOX - 1; i++)
    rPSI[i] = gsl_spline_eval (rPsispline, PSI[i], rPsiacc);
  
  // ---------------------------
  // Interpolate PSI onto r grid
  // ---------------------------
  gsl_spline_init (PSIrspline, rPSI, PSI, NRBOX);

  // ---------------------------
  // Calculate profile functions
  // ---------------------------
  double f1c   = 1./qc;
  double f3c   = - f1c * (HHfunc(2, 0) * HHfunc(2, 0) + VVfunc(2, 0) * VVfunc(2, 0));
  double p2ppc = - 2. * pc * mu;
  double g2pc  = - 2. * (f1c*f1c + p2ppc/2.);

  T  [0] = B0EXP * R0EXP;
  TTp[0] = B0EXP                      * g2pc  /(f1c + eps*eps * f3c);
  P  [0] = (B0EXP*B0EXP /mu_0)        * epsa*epsa * pc;
  Pp [0] = (B0EXP /mu_0 /R0EXP/R0EXP) * p2ppc /(f1c + eps*eps * f3c);
  Q  [0] = q2[0];
  
  for (int i = 1; i < NRBOX; i++)
    {
      double rp = rPSI[i];

      double f1  = Getf1 (rp);
      double f1p = Getf1p(rp);
      double p2  = Getp2 (rp);
      double p2p = Getp2p(rp);

      double g2  = gsl_spline_eval (g2spline, rp, g2acc);
      double f   = gsl_spline_eval (fspline,  rp, facc);
      double q   = gsl_spline_eval (q2spline, rp, q2acc);

      double g   = 1. + epsa*epsa * g2;
      double g2p = - p2p - f1*f1p/rp/rp;

      T  [i] = B0EXP * R0EXP              * g;
      TTp[i] = B0EXP                      * rp * g2p * g /f;
      P  [i] = (B0EXP*B0EXP /mu_0)        * epsa*epsa * p2;
      Pp [i] = (B0EXP /mu_0 /R0EXP/R0EXP) * rp * p2p     /f;
      Q  [i] = q;                             
    }

  // -------------------------
  // Calculate boundary points
  // -------------------------
  NPBOUND = Nw + 1;
  for (int i = 0; i < NPBOUND; i++)
    {
      RBOUND[i] = R0EXP * Rbound[i];
      ZBOUND[i] = R0EXP * Zbound[i];
    }

  // -----------------------------------------
  // Calculate bounding box and limiter points
  // -----------------------------------------
  double Rmin = 1.e6, Rmax = -1.e6, Zmin = 1.e6, Zmax = -1.e6;

  for (int i = 0; i < NPBOUND; i++)
    {
      if (RBOUND[i] < Rmin)
	Rmin = RBOUND[i];
      if (RBOUND[i] > Rmax)
	Rmax = RBOUND[i];
      if (ZBOUND[i] < Zmin)
	Zmin = ZBOUND[i];
      if (ZBOUND[i] > Zmax)
	Zmax = ZBOUND[i];
    }
  
  Rmin = Rmin - 0.05 * (Rmax - Rmin);
  Rmax = Rmax + 0.05 * (Rmax - Rmin);
  Zmin = Zmin - 0.05 * (Zmax - Zmin);
  Zmax = Zmax + 0.05 * (Zmax - Zmin);

  RBOXLFT = Rmin;
  RBOXLEN = Rmax - Rmin;
  ZOFF    = (Zmin + Zmax) /2.;
  ZBOXLEN = Zmax - Zmin;

  NLIMITER = 5;

  RLIMITER[0] = Rmin; ZLIMITER[0] = Zmin;
  RLIMITER[1] = Rmax; ZLIMITER[1] = Zmin;
  RLIMITER[2] = Rmax; ZLIMITER[2] = Zmax;
  RLIMITER[3] = Rmin; ZLIMITER[3] = Zmax;
  RLIMITER[4] = Rmin; ZLIMITER[4] = Zmin;

  // ------------------------------------
  // Set up normalized R and Z gridpoints
  // ------------------------------------
  for (int i = 0; i < NRBOX; i++)
    RGRID[i] = (RBOXLFT           + RBOXLEN * double (i) /double (NRBOX - 1)) /R0EXP;
  for (int j = 0; j < NZBOX; j++)
    ZGRID[j] = (ZOFF - ZBOXLEN/2. + ZBOXLEN * double (j) /double (NZBOX - 1)) /R0EXP;

  // -----------------------
  // Calculate PSIRZ on grid
  // -----------------------
  int cnt = 0;
  for (int i = 0; i < NRBOX; i++)
    {
      for (int j = 0; j < NZBOX; j++)
	{
	  double R = RGRID[i];
	  double Z = ZGRID[j];
	  
	  double r = sqrt ((R - 1.) * (R - 1.) + Z*Z) /epsa;
	  double w = atan2 (Z, 1. - R);
	  
	  for (int k = 0; k < 10; k++)
	    {
	      double RR = R - Getf_R (r, w);
	      double ZZ = Z - Getf_Z (r, w);
	      
	      r = sqrt ((RR - 1.) * (RR - 1.) + ZZ*ZZ) /epsa;
	      w = atan2 (ZZ, 1. - RR);
	    }
	  rRZ[cnt] = r;
	  wRZ[cnt] = w;
	  
	  if (r > rc)
	    PSIRZ[cnt] = epsa*epsa * B0EXP*R0EXP*R0EXP * GetPSIvac (rc) * r*r/rc/rc;
	  else if (r >= 1.)
	    PSIRZ[cnt] = epsa*epsa * B0EXP*R0EXP*R0EXP * GetPSIvac (r);
	  else
	    PSIRZ[cnt] = gsl_spline_eval (PSIrspline, r, PSIracc);
	  
	  cnt++;
	}
      if (i%128 == 0)
	{
	  printf ("%4d ", i); fflush (stdout);
	}
    }
  printf ("\n");

  // --------------------------------------
  // Set up unnormalized R and Z gridpoints
  // --------------------------------------
  for (int i = 0; i < NRBOX; i++)
    RGRID[i] *= R0EXP;
  for (int j = 0; j < NZBOX; j++)
    ZGRID[j] *= R0EXP;

  // -------------------------------
  // Output EFIT data to netcdf file
  // -------------------------------
  printf ("Writing EFIT data to netcdf file Outputs/Equilibrium/EFIT.nc:\n");
  
  int pint[4];

  pint[0] = NRBOX;
  pint[1] = NZBOX;
  pint[2] = NPBOUND;
  pint[3] = NLIMITER;

  double preal[15];

  preal[0]  = RBOXLEN;
  preal[1]  = ZBOXLEN;
  preal[2]  = RBOXLFT;
  preal[3]  = ZOFF;
  preal[4]  = R0EXP;
  preal[5]  = B0EXP;
  preal[6]  = RAXIS;
  preal[7]  = ZAXIS;
  preal[8]  = PSIAXIS;
  preal[9]  = PSIBOUND;
  preal[10] = CURRENT;
  preal[11] = Rmin;
  preal[12] = Rmax;
  preal[13] = Zmin;
  preal[14] = Zmax;
  
  try
    {
      NcFile dataFile ("Outputs/Equilibrium/EFIT.nc", NcFile::replace);

      NcDim i_d = dataFile.addDim ("Ni", 4);
      NcDim r_d = dataFile.addDim ("Nr", 15);
      NcDim p_d = dataFile.addDim ("Np", NRBOX);
      NcDim z_d = dataFile.addDim ("Nz", NZBOX);
      NcDim b_d = dataFile.addDim ("Nb", NPBOUND);
      NcDim l_d = dataFile.addDim ("Nl", NLIMITER);

      vector<NcDim> psi_d;
      psi_d.push_back (p_d);
      psi_d.push_back (z_d);
 
      NcVar i_x   = dataFile.addVar ("IntegerParameters", ncInt,    i_d);
      i_x.putVar (pint);
      NcVar r_x   = dataFile.addVar ("RealParameters",    ncDouble, r_d);
      r_x.putVar (preal); 
      NcVar pn_x  = dataFile.addVar ("PSI_N",             ncDouble, p_d);
      pn_x.putVar (PSIN);
      NcVar T_x   = dataFile.addVar ("T",                 ncDouble, p_d);
      T_x.putVar (T);
      NcVar P_x   = dataFile.addVar ("P",                 ncDouble, p_d);
      P_x.putVar (P);
      NcVar TTp_x = dataFile.addVar ("TTp",               ncDouble, p_d);
      TTp_x.putVar (TTp);
      NcVar Pp_x  = dataFile.addVar ("Pp",                ncDouble, p_d);
      Pp_x.putVar (Pp);
      NcVar PSI_x = dataFile.addVar ("PSI",               ncDouble, psi_d);
      PSI_x.putVar (PSIRZ);
      NcVar rx_x = dataFile.addVar  ("r",                 ncDouble, psi_d);
      rx_x.putVar (rRZ); 
      NcVar w_x = dataFile.addVar   ("w",                 ncDouble, psi_d);
      w_x.putVar (wRZ);
      NcVar Q_x   = dataFile.addVar ("Q",                 ncDouble, p_d);
      Q_x.putVar (Q);
      NcVar R_x   = dataFile.addVar ("RBOUND",            ncDouble, b_d);
      R_x.putVar (RBOUND);
      NcVar Z_x   = dataFile.addVar ("ZBOUND",            ncDouble, b_d);
      Z_x.putVar (ZBOUND);
      NcVar rr_x  = dataFile.addVar ("RLIMITER",          ncDouble, l_d);
      rr_x.putVar (RLIMITER);
      NcVar zz_x  = dataFile.addVar ("ZLIMITER",          ncDouble, l_d);
      zz_x.putVar (ZLIMITER);
      NcVar rg_x  = dataFile.addVar ("RGRID",             ncDouble, p_d);
      rg_x.putVar (RGRID);
      NcVar zg_x  = dataFile.addVar ("ZGRID",             ncDouble, z_d);
      zg_x.putVar (ZGRID);
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Equilibrium/EFIT.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  // -------------------
  // Output PSI in ascii
  // -------------------
  FILE* file = OpenFilew ("Outputs/Equilibrium/PsiSequential.txt");
  cnt = 0;
  for (int i = 0; i < NRBOX; i++)
    for (int j = 0; j < NRBOX; j++)
      {
	fprintf (file, "%4d %4d %17.9e\n", i, j, PSIRZ[cnt]);
	cnt++;
      }
  fclose (file);

  // --------
  // Clean up
  // --------
  delete[] PSI;    delete[] rPSI;   delete[] T;        delete[] TTp;      delete[] Pp;   
  delete[] RBOUND; delete[] ZBOUND; delete[] RLIMITER; delete[] ZLIMITER; delete[] PSIRZ;
  delete[] RGRID;  delete[] ZGRID;  delete[] P;        delete[] PSIN;     delete[] rRZ;
  delete[] wRZ;    delete[] Q;

  gsl_spline_free (rPsispline);    gsl_spline_free (PSIrspline);
  gsl_interp_accel_free (rPsiacc); gsl_interp_accel_free (PSIracc);
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

