// LightEquilibrium.cpp

#include "LightEquilibrium.h"

// ###########
// Constructor
// ###########
LightEquilibrium::LightEquilibrium (double _qc, double _epsa, double _pc, vector<double>& _Hna, vector<double>& _Vna)
{
  // ---------------------------
  // Set root finding parameters
  // ---------------------------
  Maxiter = 60;

  // ----------------------------------------------------
  // Set control parameters passed from class Equilibrium
  // ----------------------------------------------------
  qc   = _qc;
  epsa = _epsa;
  pc   = _pc;
  for (int i = 0; i < _Hna.size(); i++)
    {
      Hna.push_back (_Hna[i]);
      Vna.push_back (_Vna[i]);
    }
  
  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/Equilibrium.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  mu      = JSONData["mu"]  .get<double> ();
  eps     = JSONData["eps"] .get<double> ();
  Ns      = JSONData["Ns"]  .get<int>    ();
  Nr      = JSONData["Nr"]  .get<int>    ();
  acc     = JSONData["acc"] .get<double> ();
  h0      = JSONData["h0"]  .get<double> ();
  hmin    = JSONData["hmin"].get<double> ();
  hmax    = JSONData["hmax"].get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (qc < 0.)
    {
      printf ("LightEquilibrium:: Error - qc cannot be negative\n");
      exit (1);
    }
  if (pc < 0.)
    {
      printf ("LightEquilibrium:: Error - pc cannot be negative\n");
      exit (1);
    }
  if (mu < 1.)
    {
      printf ("LightEquilibrium:: Error - mu cannot be less than unity\n");
      exit (1);
    }
  if (epsa < 0.)
    {
      printf ("LightEquilibrium:: Error - epsa cannot be less than unity\n");
      exit (1);
    }
  if (eps < 0.)
    {
      printf ("LightEquilibrium:: Error - eps cannot be less than unity\n");
      exit (1);
    }
  if (Ns < 1)
    {
      printf ("LightEquilibrium:: Error - Ns cannot be less than unity\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("LightEquilibrium:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("LightEquilibrium:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("LightEquilibrium:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("LightEquilibrium:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("LightEquilibrium:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("LightEquilibrium:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (Hna.size() != Vna.size())
    {
      printf ("LightEquilibrium:: Error - Hna and Vna arrays must be the same size\n");
      exit (1);
    }
 }

// ##########
// Destructor
// ##########
LightEquilibrium::~LightEquilibrium ()
{
}

// #####################################################################
// Function to calculate nu value that gives required edge safety-factor
// #####################################################################
void LightEquilibrium::GetNu (double _qa, double& _nu)
{
  qa = _qa;

  double nustart = 0.5 * qa /qc;
  double nuend   = 2.  * qa /qc;

  _nu = RootFind (nustart, nuend);

  double qcentral, qedge, res, sa, sat;
  GetSafety (_nu, qcentral, qedge, sa, sat);
  res = fabs (qedge - qa);
  
  printf ("\nClass LIGHTEQUILIBRIUM::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("nu     = %10.3e q_central     = %10.3e q_edge = %10.3e res = %10.3e\n",
	  _nu, qcentral, qedge, res);
  printf ("s_edge = %10.3e s_edge_target = %10.3e                     res = %10.3e\n",
	  sa, sat, fabs (sa - sat));

  if (res > 1.e-3)
    {
      printf ("LightEquilibrium:: Error - search residual too large: %10.3e\n", res);
      exit (1);
    }
}

// ############################################
// Function to return central and edge q-values
// ############################################
void LightEquilibrium::GetSafety (double _nu, double& qcentral, double& qedge, double& sa, double& sat)
{
  nu = _nu;
  
  // ...............
  // Allocate memory
  // ...............
  rr = new double[Nr+1];
  f1 = new double[Nr+1];
  f3 = new double[Nr+1];
  g2 = new double[Nr+1];
  q0 = new double[Nr+1];
  q2 = new double[Nr+1];

  HHfunc.resize (Ns+1, Nr+1);
  VVfunc.resize (Ns+1, Nr+1);
  HPfunc.resize (Ns+1, Nr+1);
  VPfunc.resize (Ns+1, Nr+1);

  g2spline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

  g2acc = gsl_interp_accel_alloc ();

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

  // ..................
  // Set up radial grid
  // ..................
  for (int i = 0; i <= Nr; i++)
    {
      rr[i] = double (i) /double (Nr);
      f1[i] = Getf1  (rr[i]);
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
  
  // .........................
  // Rescale shaping functions
  // .........................
  int    nn   = 1;
  double zero = 0.;
 
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

  // .....................
  // Integrate f3 equation
  // .....................
  double* y1    = new double[1];
  double* dy1dr = new double[1];
  double* err1  = new double[1];
  rhs_chooser   = 1;

  f1c        = 1. /qc;
  double H2c = HPfunc(2, 0);
  double V2c = VPfunc(2, 0);

  h       = h0;
  count   = 0;
  f3[0]   = 0.;
  q0[0]   = qc;

  r     = eps;
  y1[0] = - f1c * (H2c*H2c + V2c*V2c) * r*r;
  
  for (int i = 1; i <= Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (1, r, y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i] - h);
      CashKarp45Fixed (1, r, y1, err1, rr[i] - r);
      CashKarp45Rhs (rr[i], y1, dy1dr);
      
      f3[i] = y1[0];
      q0[i] = rr[i]*rr[i] /f1[i];
      q2[i] = rr[i]*rr[i] * (1. + epsa*epsa*g2[i]) * exp(- epsa*epsa * f3[i]/f1[i]) /f1[i];
    }
  q2[0] = qc * (1. + epsa*epsa * (H2c*H2c + V2c*V2c));

  double ff1 = f1[Nr];
  double ff3 = f3[Nr];
  double f1p = Getf1p (1.);
  double p2p = Getp2p (1.);
  double gg2 = g2[Nr];
  double g2p = - p2p - ff1*f1p;
  double f3p = dy1dr[0];
  
  sa = 2. - epsa*epsa * f3p/ff1;

  delete[] y1; delete[] dy1dr; delete[] err1;

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
  sat = 2. + epsa*epsa * sum;

  // ........
  // Clean up
  // ........
  delete[] rr;  delete[] f1; delete[] f3; delete[] g2; delete[] q0;  delete[] q2; 
  
  gsl_spline_free (g2spline);
 
  gsl_interp_accel_free (g2acc);
 
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

  qcentral = q2[0];
  qedge    = q2[Nr];
}

// ########################
// Function to return f1(r)
// ########################
double LightEquilibrium::Getf1 (double r)
{
  if (r < 0.1)
    return r*r * (1. - (nu-1.)*r*r/2. + (nu-1.)*(nu-2.)*r*r*r*r/6. - (nu-1.)*(nu-2.)*(nu-3.)*r*r*r*r*r*r/24.)/qc;
  else
    return (1. - pow (1. - r*r, nu)) /nu/qc;
}

// #########################
// Function to return f1'(r)
// #########################
double LightEquilibrium::Getf1p (double r) 
{
  if (r < 0.1)
    return 2.*r * (1. - (nu-1.)*r*r + (nu-1.)*(nu-2.)*r*r*r*r/2. - (nu-1.)*(nu-2.)*(nu-3.)*r*r*r*r*r*r/6.)/qc;
  else
    return 2.*r * pow (1. - r*r, nu-1.) /qc;
}

// #########################
// Function to return p2'(r)
// #########################
double LightEquilibrium::Getp2p (double r)
{
  return - 2. * pc * mu * r * pow (1. - r*r, mu-1.);
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void LightEquilibrium::CashKarp45Rhs (double r, double* y, double* dydr)
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

      delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;
    }
 }

// ##############################################
// Target function for 1-dimensional root finding
// ##############################################
double LightEquilibrium::RootFindF (double x)
{
  double qcentral, qedge, sa, sat;

  GetSafety (x, qcentral, qedge, sa, sat);

  return qedge - qa;
}

