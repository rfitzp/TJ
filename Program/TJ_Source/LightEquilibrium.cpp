// LightEquilibrium.cpp

#include "LightEquilibrium.h"

// ###########
// Constructor
// ###########
LightEquilibrium::LightEquilibrium ()
{
  // -----------------------------------
  // Set adaptive integration parameters
  // -----------------------------------
  maxrept = 50;
  flag    = 2;

  // ---------------------------
  // Set root finding parameters
  // ---------------------------
  nint    = 10;
  Eta     = 1.e-12;
  Maxiter = 60;

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

  qc   = JSONData["Equilibrium_control"]["qc"]  .get<double> ();
  pc   = JSONData["Equilibrium_control"]["pc"]  .get<double> ();
  mu   = JSONData["Equilibrium_control"]["mu"]  .get<double> ();
  epsa = JSONData["Equilibrium_control"]["epsa"].get<double> ();
  eps  = JSONData["Equilibrium_control"]["eps"] .get<double> ();
  Ns   = JSONData["Equilibrium_control"]["Ns"]  .get<int> ();
  Nr   = JSONData["Equilibrium_control"]["Nr"]  .get<int> ();
  acc  = JSONData["Equilibrium_control"]["acc"] .get<double> ();
  h0   = JSONData["Equilibrium_control"]["h0"]  .get<double> ();
  hmin = JSONData["Equilibrium_control"]["hmin"].get<double> ();
  hmax = JSONData["Equilibrium_control"]["hmax"].get<double> ();

  for (const auto& number : JSONData["Equilibrium_control"]["Hna"])
    {
      Hna.push_back (number.get<double> ());
    }
  for (const auto& number : JSONData["Equilibrium_control"]["Vna"])
    {
      Vna.push_back (number.get<double> ());
    }

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
      printf ("LightEquilibrium:: Error - Hna and Van arrays must be the same size\n");
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
double LightEquilibrium::GetNu (double _qa)
{
  qa = _qa;

  double nustart = 0.5 * qa /qc;
  double nuend   = 2.  * qa /qc;

  double _nu = RootFind (nustart, nuend);

  double qcentral, qedge, res;
  GetSafety (_nu, qcentral, qedge);
  res = fabs (qedge - qa);

  printf ("\nClass LIGHTEQUILIBRIUM::\n");
  printf ("nu = %10.3e q_central = %10.3e q_edge = %10.3e res = %10.3e\n", _nu, qcentral, qedge, res);

  if (res > 1.e-3)
    {
      printf ("LightEquilibrium:: Error - search residual too large: %10.3e\n", res);
      exit (1);
    }

  return _nu;
}

// ############################################
// Function to return central and edge q-values
// ############################################
void LightEquilibrium::GetSafety (double _nu, double& qcentral, double& qedge)
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
      rr [i] = double (i) /double (Nr);
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
      Rhs (rr[i], y1, dy1dr);
      
      f3[i] = y1[0];
      q0[i] = rr[i]*rr[i] /f1[i];
      q2[i] = rr[i]*rr[i] * (1. + epsa*epsa*(g2[i] - f3[i]/f1[i])) /f1[i];

    }
  q2[0] = qc * (1. + epsa*epsa * (H2c*H2c + V2c*V2c));

  delete[] y1; delete[] dy1dr; delete[] err1;

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
    return 2. * r * pow (1. - r*r, nu-1.) /qc;
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
void LightEquilibrium::Rhs (double r, double* y, double* dydr)
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
void LightEquilibrium::CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
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
void LightEquilibrium::CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h)
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

// ################################
// Target function for zero finding
// ################################
double LightEquilibrium::Feval (double x)
{
 double qcentral, qedge;

 GetSafety (x, qcentral, qedge);

 return qedge - qa;
}

// ###################################################################
// Routine to find approximate root of F(x) = 0 using Ridder's method 
// Search takes place in interval (x1, x2)
// Interval is chopped into nint equal segments
// ###################################################################
double LightEquilibrium::RootFind (double x1, double x2)
{
  double F1, F2 = 0., root;

  // Chop search interval into nint segments  
  for (int i = 0; i < nint; i++)
    {
      double x1_seg = x1 + (x2 - x1) * double (i)     /double (nint);
      double x2_seg = x1 + (x2 - x1) * double (i + 1) /double (nint);
      
      if (i == 0) 
	F1 = Feval (x1_seg);
      else 
	F1 = F2;
      F2 = Feval (x2_seg);
      //printf ("%e %e %e %e\n", x1_seg, F1, x2_seg, F2);
      
      // Call Ridder's method for segment containing zero
      if (F1 * F2 < 0.)
	{
	  Ridder (x1_seg, x2_seg, F1, F2, root);
	  break;
	}
    }

  return root;
}

// ############################################
// Ridder's method for finding root of F(x) = 0
// ############################################
void LightEquilibrium::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) /2.;    
      double F3 = Feval (x3);
      
      // Iterate x using Ridder's method 
      xold = x;           
      x = x3 - (x3 - x1) * (F2 - F1) * F3 /
	(sqrt (F3 * F3 - F1 * F2) * fabs (F2 - F1));
      Fx = Feval (x);
       
      // Make new value of x upper/lower bound of refined search interval, as appropriate 
      if (Fx * F1 < 0.) 
	{  
	  x2 = x;           
	  F2 = Fx; 
	}
      else 
	{
	  x1 = x;
	  F1 = Fx; 
	}
      //printf ("%d %e %e\n", iter, x, Fx);
      iter++;
    } 
  // Iterate until absolute change in x falls below Eta
  while (fabs (x - xold) > Eta && fabs(Fx) > Eta && iter < Maxiter); 
}

// ##########################
// Function to read JSON file
// ##########################
json LightEquilibrium::ReadJSONFile (const string& filename)
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

