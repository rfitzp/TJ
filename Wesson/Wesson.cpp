// Wesson.cpp

#include "Wesson.h"

// ###########
// Constructor
// ###########
Wesson::Wesson ()
{
  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "Wesson.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  mpol = JSONData["mpol"].get<double> ();
  ntor = JSONData["ntor"].get<double> ();
  q0   = JSONData["q0"]  .get<double> ();
  qa   = JSONData["qa"]  .get<double> ();
  rw   = JSONData["rw"]  .get<double> ();
  abar = JSONData["abar"].get<double> ();
  qafc = JSONData["qafc"].get<double> ();
  flag = JSONData["flag"].get<int>    ();

  eps   = JSONData["eps"].get<double> ();
  del   = JSONData["del"].get<double> ();
 
  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  Eta     = JSONData["Eta"]    .get<double> ();
  Maxiter = JSONData["Maxiter"].get<int>    ();
  Nint    = JSONData["Nint"]   .get<int>    ();

  printf ("\nProgram WESSON:\n");
  printf ("mpol = %2d  ntor = %2d  q0 = %10.3e  qa = %10.3e  rw = %10.3e  abar = %10.3e  flag = %1d\n\n",
	  int(mpol), int(ntor), q0, qa, rw, abar, flag);
}

// #########################
// Function to solve problem
// #########################
void Wesson::Solve ()
{
  int    N;
  double q0start, q0end, q0crit, qastart, qaend;
  
  if (flag == 0)
    {
      // #####################
      // Stability limit scans
      // #####################
      
      N       = 1000;
      qastart =        mpol /ntor + 1.e-3;
      qaend   = qafc * mpol /ntor - 1.e-3;
      
      FILE* file1 = fopen ("Output/Stability.out", "w");
      for (int i = 0; i < N; i++)
	{
	  qa = qastart + (qaend - qastart) * double (i) /double (N-1);
	  
	  q0start = 0.5;
	  q0end   = qa /2./abar/abar;
	  if (q0end > mpol/ntor - 1.e-3)
	    q0end = mpol /ntor - 1.e-3;
	  
	  q0crit = RootFind (q0start, q0end);
	  
	  double nucrit = qa /q0crit /abar/abar - 1.;
	  double li     = Getli ();
	  double li1    = Getli1 ();

	  if (q0crit < 0.)
	    break;
	  
	  printf ("qa = %10.3e q0_start = %10.3e q0_end = %10.3e q0_crit = %10.3e nu_crit_0 = %10.3e li = %10.3e\n",
		  qa, q0start, q0end, q0crit, nucrit, li);
	  fprintf (file1, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", qa, nucrit, qa/(mpol/ntor)/abar/abar - 1., li, li1);
	  fflush (file1);
	}
      
      fclose (file1);
    }
  else
    {
      // ####################################
      // 2/1 mode Delta' scans at constant q0
      // ####################################

      N       = 1000;
      q0      = 0.8;
      mpol    = 2.;
      ntor    = 1.;
      abar    = 1.;
      qastart = 2.0001;
      qaend   = 20.;
      
      // %%%%%%%%
      // rw = 1.2
      // %%%%%%%%
      FILE* file = fopen ("Output/m2n1r12.out", "w");
      
      rw = 1.2;
      
      for (int i = 0; i < N; i++)
	{
	  qa = qastart + (double (i) / double (N)) * (qaend - qastart);
	  
	  double Deltas = GetDelta ();
	  
	  double rs     = Findrs (mpol /ntor);
	  double alphas = Getalpha (rs);
	  double betas  = Getbeta (rs);
	  double Wsat;
	  if (Deltas > 0.)
	    Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
	  else
	    Wsat = 0.;
	  double wsat = rs / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas) /4.;
	  double ss   = Gets(rs);
	  
	  printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e wsat = %11.4e  ss = %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	}
      
      fclose (file);
      
      // %%%%%%%%
      // rw = 1.1
      // %%%%%%%%
      file = fopen ("Output/m2n1r11.out", "w");
      
      rw = 1.1;
      
      for (int i = 0; i < N; i++)
	{
	  qa = qastart + (double (i) / double (N)) * (qaend - qastart);
	  
	  double Deltas = GetDelta ();
	  
	  double rs     = Findrs (mpol /ntor);
	  double alphas = Getalpha (rs);
	  double betas  = Getbeta (rs);
	  double Wsat;
	  if (Deltas > 0.)
	    Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
	  else
	    Wsat = 0.;
	  double wsat = rs / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas) /4.;
	  double ss   = Gets(rs);
	  
	  printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e wsat = %11.4e  ss = %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	}
      
      fclose (file);
      
      // %%%%%%%%
      // rw = 1.0
      // %%%%%%%%
      file = fopen ("Output/m2n1r10.out", "w");
      
      rw = 1.001;
      
      for (int i = 0; i < N; i++)
	{
	  qa = qastart + (double (i) / double (N)) * (qaend - qastart);
	  
	  double Deltas = GetDelta ();
	  
	  double rs     = Findrs (mpol /ntor);
	  double alphas = Getalpha (rs);
	  double betas  = Getbeta (rs);
	  double Wsat;
	  if (Deltas > 0.)
	    Wsat = rs * Deltas / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas);
	  else
	    Wsat = 0.;
	  double wsat = rs / (0.8 * alphas*alphas - 0.27 * betas  - 0.09 * alphas) /4.;
	  double ss   = Gets(rs);
	  
	  printf ("qa = %11.4e  rs = %11.4e  Delta = %11.4e  Wsat = %11.4e wsat = %11.4e  ss = %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", qa, rs, Deltas, Wsat, wsat, ss);
	}
      
      fclose (file);
    }
}

// #############################################
// Function to calculate tearing stability index
// #############################################
double Wesson::GetDelta ()
{
  // Find radius of resonant surface
  double qs = mpol /ntor;
  double rs = Findrs (qs);

  // Launch solution from magnetic axis and integrate to resonant surface
  double          r, h, t_err;
  int             rept; count = 0;
  Array<double,1> y(2);

  r    = eps;
  y(0) = pow (r, mpol);
  y(1) = mpol * pow (r, mpol - 1.);
  h    = h0;
  flg  = 0;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r < rs - del);
  RK4Fixed (r, y, rs - del - r);

  double alphas = Getalpha (r);
  double Deltam = r * y(1) /y(0) - alphas * (1. + log(del));

  //printf ("rm = %11.4e  rs - rm = %11.4e  Deltam = %11.4e\n", r, rs - r, Deltam);

  // Launch solution from plasma boundary and integrate to resonant surface

  double fac = pow (1./rw, 2.*mpol);
  r    = 1.;
  y(0) = 1.; 
  y(1) = - mpol * (1. + fac) /(1. - fac);
  h    = - h0;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r > rs + del);
  RK4Fixed (r, y, rs + del - r);

  alphas = Getalpha (r);
  double Deltap = r * y(1) /y(0) - alphas * (1. + log(del));

  //printf ("rm = %11.4e  rs - rm = %11.4e  Deltap = %11.4e\n", r, rs - r, Deltap);

  return Deltap - Deltam;
}

// ########################################
// Function to find resonant surface radius
// ########################################
double Wesson::Findrs (double qs)
{
  if (qs <= q0 || qs >= qa)
    {
      printf ("Wesson:Findrs - Error: no resonant surface in plasma\n");
      exit (1); 
    }

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
	  if (i == 0)
	    {
	      double nu = Getnu ();

	      return abar * sqrt ((qs/q0 - 1.) * (nu/2.));
	    }
	  else
	    {
	      double rn = r1, qn;
	      int count = 0;
	      do
		{
		  qn        = Getq (rn);
		  double sn = Gets (rn);
		  
		  //printf ("rn = %11.4e  qn = %11.4e  residual = %11.4e\n", rn, qn, fabs (qn - qs));
		  
		  rn += rn * (qs /qn - 1.) /sn;
		  count++;
		}
	      while (fabs (qn - qs) > 1.e-15 && count < 100);
	      
	      if (isnan (rn))
		{
		  printf ("Wesson:Findrs - Error: rn is NaN\n");
		  exit (1);
		}
	      else
		return rn;
	    }
	}
    }
}

// ########################################
// Function to return value of parameter nu
// ########################################
double Wesson::Getnu ()
{
  return qa /q0 /abar/abar - 1.;
}

// ##########################################
// Function to return value of plasma current
// ##########################################
double Wesson::GetJ (double r)
{
  if (r < abar)
    {
      double nu = Getnu ();
      double r2 = r*r /abar/abar;
      
      return 2. * pow (1. - r2, nu) / q0;
    }
  else
    return 0.;
}

// ##############################################################
// Function to return value of first derivative of plasma current
// ##############################################################
double Wesson::GetJp (double r)
{
  if (r < abar)
    {
      double nu = Getnu ();
      double r2 = r*r /abar/abar;
      
      return - 4. * nu * r * pow (1. - r2, nu - 1.) /q0 /abar/abar;
    }
  else
    return 0.;
}

// ###############################################################
// Function to return value of second derivative of plasma current
// ###############################################################
double Wesson::GetJpp (double r)
{
  if (r < abar)
    {
      double nu = Getnu ();
      double r2 = r*r /abar/abar;

      return - 4. * nu * pow (1. - r2, nu - 2.) * (1. - (2.*nu - 1.) * r2) /q0 /abar/abar;
    }
  else
    return 0.;
}

// #########################################
// Function to return value of safety-factor
// #########################################
double Wesson::Getq (double r)
{
  double nu = Getnu ();
  double r2 = r*r /abar/abar;

  if (r < eps)
    return q0 * (1. + 0.5 * nu * r2);
  else if (r < abar)
    return qa * abar*abar * r2 /(1. - pow (1. - r2, nu + 1.));
  else
    return qa * abar*abar * r2;
}

// ##########################################
// Function to return value of magnetic shear
// ##########################################
double Wesson::Gets (double r)
{
  double J = GetJ (r);
  double q = Getq (r);

  return 2. - q * J; 
}

// ###########################################
// Function to return value of parameter alpha
// ###########################################
double Wesson::Getalpha (double r)
{
  double Jp = GetJp (r);
  double q  = Getq  (r);
  double s  = Gets  (r);

  return - q * r * Jp /s; 
}

// ##########################################
// Function to return value of parameter beta
// ##########################################
double Wesson::Getbeta (double r)
{
  double Jpp = GetJpp (r);
  double q   = Getq   (r);
  double s   = Gets   (r);

  return - q * r*r * Jpp /s; 
}

// ######################################################
// Function to return value of normalized self-inductance
// ######################################################
double Wesson::Getli ()
{
  double          r, h, t_err;
  int             rept; count = 0;
  Array<double,1> y(1);

  r    = eps;
  y(0) = 0.;
  h    = h0;
  flg  = 1;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r < 1.);
  RK4Fixed (r, y, 1. - r);

  return y(0) + log (1. /abar/abar);
}

// ######################################################
// Function to return value of normalized self-inductance
// ######################################################
double Wesson::Getli1 ()
{
  double          r, h, t_err;
  int             rept; count = 0;
  Array<double,1> y(1);

  r    = eps;
  y(0) = 0.;
  h    = h0;
  flg  = 2;

  do
    {
      RK4Adaptive (r, y, h, t_err, acc, 2., rept, 20, hmin, 1.e-2, 2, 0, NULL);
    }
  while (r < 1.);
  RK4Fixed (r, y, 1. - r);

  return y(0) + log (1. /abar/abar);
}

// ##############################################################
// Function to evaluate right-hand side of differential equations
// ##############################################################
void Wesson::Rhs (double r, Array<double,1>& y, Array<double,1>& dydr)
{
  if (flg == 0)
    {
      dydr(0) = y(1);
      dydr(1) = - y(1) /r + mpol*mpol * y(0) /r/r + GetJp (r) * y(0) /r /(1./Getq (r) - ntor /mpol);
    }
  else if (flg == 1)
    {
      double nu = Getnu ();

      if (r < 1.)
	dydr(0) = (1. - pow (1. - r, nu + 1.)) * (1. - pow (1. - r, nu + 1.)) /r;
      else
	dydr(0) = 1. /r;
    }
   else
    {
      double nu = qa /(mpol/ntor) /abar/abar - 1.;

      if (r < 1.)
	dydr(0) = (1. - pow (1. - r, nu + 1.)) * (1. - pow (1. - r, nu + 1.)) /r;
      else
	dydr(0) = 1. /r;
    }
}
 
// ######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive fourth-order Runge-Kutta scheme
//
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... step-length cannot change by more than this factor from
//                  step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step-length
//     h_max   ... maximum allowable step-length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step whilst attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
//  If step-length falls below h_min then routine aborts
// ######################################################################
void Wesson::RK4Adaptive (double& x, Array<double,1>& y, double& h, 
			     double& t_err, double acc, double S, int& rept,
			     int maxrept, double h_min, double h_max, int flag, 
			     int diag, FILE* file)
{
  int neqns = y.extent(0);
  Array<double,1> y0(neqns), y1(neqns);
  double hin = h;

  // Save initial data
  double x0 = x;
  y0 = y;

  // Take full step 
  RK4Fixed (x, y, h);

  // Save data
  y1 = y;

  // Restore initial data 
  x = x0;
  y = y0;

  // Take two half-steps 
  RK4Fixed (x, y, h/2.);
  RK4Fixed (x, y, h/2.);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (y(i) - y1(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs ((y(i) - y1(i)) / y(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs ((y(i) - y1(i)) / y(i));
          err2  =  fabs (y(i) - y1(i));
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15) t_err = 1.e-15;

  // Calculate new step-length 
  double h_est = h * pow (fabs (acc / t_err), 0.2);

  // Prevent step-length from changing by more than factor S
  if (h_est / h > S)
    h *= S;
  else if (h_est / h < 1. / S)
    h /= S;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Abort if step-length falls below h_min
  if (fabs(h) < h_min)
    { 
      //printf ("Wesson::RK4Adpative: Warning - |h| < hmin at x = %11.4e\n", x);
      //exit (1);
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %11.4e hin = %11.4e err = %11.4e acc = %11.4e hest = %11.4e hout = %11.4e count = %3d\n", 
	     x, hin, t_err, acc, h_est, h, count);

  // If truncation error acceptable take step 
  if ((t_err <= acc) || (count >= maxrept))
    {  
      rept  = count;
      count = 0;
    }
  // If truncation error unacceptable repeat step 
  else 
    {
      count++;
      x = x0;
      y = y0;
      RK4Adaptive (x, y, h, t_err, acc, S, rept, 
		   maxrept, h_min, h_max, flag, diag, file);
    }
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length fourth-order Runge-Kutta scheme.
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
// #####################################################################
void Wesson::RK4Fixed (double& x, Array<double,1>& y, double h)
{
  int neqns = y.extent(0);
  Array<double,1> dydx(neqns), k1(neqns), k2(neqns), k3(neqns);
  Array<double,1> k4(neqns), f(neqns);

  // Wessonth intermediate step 
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1(i) = h * dydx(i);
      f(i)  = y(i) + k1(i) / 2.;
    }

  // First intermediate step 
  Rhs (x + h / 2., f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2(i) = h * dydx(i);
      f(i)  = y(i) + k2(i) / 2.;
    }

  // Second intermediate step 
  Rhs (x + h / 2., f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3(i) = h * dydx(i);
      f(i)  = y(i) + k3(i);
    }

  // Third intermediate step 
  Rhs (x + h, f, dydx);
  for (int i = 0; i < neqns; i++)
    k4(i) = h * dydx(i);

  // Actual step 
  for (int i = 0; i < neqns; i++)
    y(i) += k1(i) / 6. + k2(i) / 3. + k3(i) / 3. + k4(i) / 6.;
  x += h;
}

// ##################################################
// Target function for one-dimensional root finding.
// Needs to be overridden in inheriting class.
// ##################################################
double Wesson::RootFindF (double x)
{
  q0 = x;
  
  double Delta = GetDelta ();

  //printf ("q0 = %10.3e nu = %10.3e Del = %10.4e\n", x, Getnu (), Del);

  return Delta;
}

// ##################################################################
// Routine to find approximate root of F(x) = 0 using Ridder's method 
//
// Search takes place in interval (x1, x2)
// Interval is chopped into Nint equal segments
//
//  Eta     ... Minimum magnitude of F at root F(x) = 0
//  Maxiter ... Maximum number of iterations
// 
// ##################################################################
double Wesson::RootFind (double x1, double x2)
{
  double F1, F2 = 0., root = -1.e15;

  // Chop search interval into Nint segments  
  for (int i = 0; i < Nint; i++)
    {
      double x1_seg = x1 + (x2 - x1) * double (i)     /double (Nint);
      double x2_seg = x1 + (x2 - x1) * double (i + 1) /double (Nint);
      
      if (i == 0) 
	F1 = RootFindF (x1_seg);
      else 
	F1 = F2;
      F2 = RootFindF (x2_seg);
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
void Wesson::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) /2.;    
      double F3 = RootFindF (x3);
      
      // Iterate x using Ridder's method 
      xold = x;           
      x = x3 - (x3 - x1) * (F2 - F1) * F3 /
	(sqrt (F3 * F3 - F1 * F2) * fabs (F2 - F1));
      Fx = RootFindF (x);
       
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
  while (fabs (x - xold) > Eta && fabs (Fx) > Eta && iter < Maxiter); 
}

// #####################
// Function to open file
// #####################
FILE* Wesson::OpenFile (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Wesson::OpenFile: Error opening data-file\n");
      exit (1);
    }
  return file;
}
