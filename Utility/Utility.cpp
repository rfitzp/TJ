// Utility.cpp

#include "Utility.h"

// ###########
// Constructor
// ###########
Utility::Utility ()
{
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

  // -----------------------------------------------------------------------
  // Set default values of Cash-Karp RK4/RK5 adaptive integration parameters
  // -----------------------------------------------------------------------
  acc     = 1.e-12;
  h0      = 1.e-2;
  hmin    = 1.e-10;
  hmax    = 1.e-2;
  maxrept = 50;
  flag    = 2;

  // -------------------------------------------------------------
  // Set default values of one-dimensional root finding parameters
  // -------------------------------------------------------------
  Nint    = 10;
  Eta     = 1.e-12;
  Maxiter = 50;

  // -------------------------------------------------------------
  // Set default values of two-dimensional root finding parameters
  // -------------------------------------------------------------
  dS      = 1.e-6;
  Smax    = 0.1;
  Smin    = 1.e-6;
  Eps     = 1.e-12;
  alpha   = 1.e-4;
  MaxIter = 50;
};

// ##########
// Destructor
// ##########
Utility::~Utility ()
{
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// for Cash-Karp RK4/RK5 adaptive integration routine.
// Needs to be overridden in inheriting class.
// ###############################################################
void Utility::CashKarp45Rhs (double x, double* y, double* dydx)
{
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
void Utility::CashKarp45Adaptive (int neqns, double& x, double* y, double& h, 
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

// ######################################################################
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
void Utility::CashKarp45Fixed (int neqns, double& x, double* y, double* err, double h)
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
  CashKarp45Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  CashKarp45Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  CashKarp45Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  CashKarp45Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  CashKarp45Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  CashKarp45Rhs (x + aa6 * h, f, dydx);
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

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// for Cash-Karp RK4/RK5 adaptive integration routine.
// Needs to be overridden in inheriting class.
// ###############################################################
void Utility::CashKarp45Rhs (double x, complex<double>* y, complex<double>* dydx)
{
}

// #######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step-length Cash-Karp fourth-order/fifth-order
//  Runge-Kutta scheme
//
//     neqns   ... number of coupled equations
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
//     flag = 0 ... error is absolute
//     flag = 1 ... error is relative
//     flag = 2 ... error is mixed
//
// #######################################################################
void Utility::CashKarp45Adaptive (int neqns, double& x, complex<double>* y, double& h, 
				  double& t_err, double acc, double S, double T, int& rept,
				  int maxrept, double h_min, double h_max, int flag, 
				  int diag, FILE* file)
{
  complex<double>* y0  = new complex<double>[neqns];
  complex<double>* Err = new complex<double>[neqns];
  double           hin = h;

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
          err   = abs (Err[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = abs (Err[i] /y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = abs (Err[i] /y[i]);
	  err2  = abs (Err[i]);
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step-length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc /t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc /t_err), 0.25);

  // Prevent step-length from changing by more than factor T
  if (h_est /h > T)
    h *= T;
  else if (h_est /h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h /fabs(h) : h;

  // Prevent step-length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = + h_min;
      else
	h = - h_min;
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

  delete[] y0; delete[] Err;
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length Cash-Karp fourth-order/fifth-order
// Runge-Kutta scheme
//
//     neqns ... number of coupled equations
//     x     ... independent variable
//     y     ... array of dependent variables 
//     err   ... array of errors
//     h     ... step-length
//     
// #####################################################################
void Utility::CashKarp45Fixed (int neqns, double& x, complex<double>* y, complex<double>* err, double h)
{
  complex<double>* dydx = new complex<double>[neqns];
  complex<double>* k1   = new complex<double>[neqns];
  complex<double>* k2   = new complex<double>[neqns];
  complex<double>* k3   = new complex<double>[neqns];
  complex<double>* k4   = new complex<double>[neqns];
  complex<double>* k5   = new complex<double>[neqns];
  complex<double>* k6   = new complex<double>[neqns];
  complex<double>* f    = new complex<double>[neqns];

  // First stage
  CashKarp45Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  CashKarp45Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  CashKarp45Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  CashKarp45Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  CashKarp45Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  CashKarp45Rhs (x + aa6 * h, f, dydx);
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

  delete[] dydx; delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6;
  delete[] f;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// for Cash-Karp RK4/RK5 adaptive integration routine.
// Needs to be overridden in inheriting class.
// ###############################################################
void Utility::CashKarp45Rhs1 (double x, complex<double>* y, complex<double>* dydx)
{
}

// #######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step-length Cash-Karp fourth-order/fifth-order
//  Runge-Kutta scheme
//
//     neqns   ... number of coupled equations
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
//     flag = 0 ... error is absolute
//     flag = 1 ... error is relative
//     flag = 2 ... error is mixed
//
// #######################################################################
void Utility::CashKarp45Adaptive1 (int neqns, double& x, complex<double>* y, double& h, 
				   double& t_err, double acc, double S, double T, int& rept,
				   int maxrept, double h_min, double h_max, int flag, 
				   int diag, FILE* file)
{
  complex<double>* y0  = new complex<double>[neqns];
  complex<double>* Err = new complex<double>[neqns];
  double           hin = h;

  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0[i] = y[i];

  // Take Cash-Karp RK4/RK5 step 
  CashKarp45Fixed1 (neqns, x, y, Err, h);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = abs (Err[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = abs (Err[i] /y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = abs (Err[i] /y[i]);
	  err2  = abs (Err[i]);
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step-length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc /t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc /t_err), 0.25);

  // Prevent step-length from changing by more than factor T
  if (h_est /h > T)
    h *= T;
  else if (h_est /h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h /fabs(h) : h;

  // Prevent step-length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = + h_min;
      else
	h = - h_min;
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
      CashKarp45Adaptive1 (neqns, x, y, h, t_err, acc, S, T, rept, 
			   maxrept, h_min, h_max, flag, diag, file);
    }

  delete[] y0; delete[] Err;
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length Cash-Karp fourth-order/fifth-order
// Runge-Kutta scheme
//
//     neqns ... number of coupled equations
//     x     ... independent variable
//     y     ... array of dependent variables 
//     err   ... array of errors
//     h     ... step-length
//     
// #####################################################################
void Utility::CashKarp45Fixed1 (int neqns, double& x, complex<double>* y, complex<double>* err, double h)
{
  complex<double>* dydx = new complex<double>[neqns];
  complex<double>* k1   = new complex<double>[neqns];
  complex<double>* k2   = new complex<double>[neqns];
  complex<double>* k3   = new complex<double>[neqns];
  complex<double>* k4   = new complex<double>[neqns];
  complex<double>* k5   = new complex<double>[neqns];
  complex<double>* k6   = new complex<double>[neqns];
  complex<double>* f    = new complex<double>[neqns];

  // First stage
  CashKarp45Rhs1 (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  CashKarp45Rhs1 (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  CashKarp45Rhs1 (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  CashKarp45Rhs1 (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  CashKarp45Rhs1 (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  CashKarp45Rhs1 (x + aa6 * h, f, dydx);
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

  delete[] dydx; delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6;
  delete[] f;
}

// ##################################################
// Target function for one-dimensional root finding.
// Needs to be overridden in inheriting class.
// ##################################################
double Utility::RootFindF (double x)
{
  return 0.;
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
double Utility::RootFind (double x1, double x2)
{
  double F1, F2 = 0., root;

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
void Utility::Ridder (double x1, double x2, double F1, double F2, double& x)
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
  while (fabs (x - xold) > Eta && fabs(Fx) > Eta && iter < Maxiter); 
}

// ############################################################################
// Target functions for two-dimensional root finding via Newton-Raphson method.
// Needs to be overridden in inheriting class.
// ############################################################################
void Utility::NewtonFunction (double x1, double x2, double& F1, double& F2)
{
}

// #################################################################################
// Function to find root of F1(x1,x2) = F2(x1,x2) = 0 via Newton-Raphson method
//
//   dS      ... Step-size for calculation of Jacobian 
//   Smax    ... Maximum step-size
//   Smin    ... Minimum step-size 
//   Eps     ... Minimum magnitude of |F1^2+F2^2|^1/2 at root F1(x1,x2) = F2(x1,x2) = 0
//   MaxIter ... Maximum number of iterations
//   alpha   ... Ensures sufficient decrease in function value
//
// #################################################################################
void Utility::NewtonRoot (double& x1, double& x2, double& Residual, int verbose)
{
  double F1, F2, J11, J12, J21, J22, det, iJ11, iJ12, iJ21, iJ22, dx1, dx2, dx, f, g1, g2, lambda;
  int    iter;

  iter = 0;
  do
    {
      NewtonFunction (x1, x2, F1, F2);

      NewtonJacobian (x1, x2, J11, J12, J21, J22);
      
      det = J11 * J22 - J12 * J21;
      
      iJ11 =   J22 /det;
      iJ12 = - J12 /det;
      iJ21 = - J21 /det;
      iJ22 =   J11 /det;
      
      dx1 = - (iJ11 * F1 + iJ12 * F2);
      dx2 = - (iJ21 * F1 + iJ22 * F2);

      dx = sqrt (dx1*dx1 + dx2*dx2);

      f = 0.5 * (F1*F1 + F2*F2);

      g1 = F1*J11 + F2*J21;
      g2 = F1*J12 + F2*J22;

      NewtonBackTrack (x1, x2, dx1, dx2, dx, f, g1, g2, lambda);
 
      NewtonFunction (x1, x2, F1, F2);
	
      Residual = sqrt (F1*F1 + F2*F2);

      if (verbose)
	printf ("x = (%10.3e, %10.3e) F = (%10.3e, %10.3e) Residual = %10.3e lambda = %10.3e dx = %10.3e\n",
		x1, x2, F1, F2, Residual, lambda, dx);
    
      iter++;
    }
  while (Residual > Eps && dx > Smin && iter < MaxIter);
}

// #################################################################################################
// Function to backtrack along Newton step in order to minimize f = (F1*F1 + F2*F2) /2
//
// Press, Teukolsky, Vetterling, and Flannery, Numerical Recipies in C (Cambridge, 1992), Sect. 9.7
//
// ##################################################################################################
void Utility::NewtonBackTrack (double& x1, double& x2, double dx1, double dx2, double& dx, double f, double g1, double g2, double& lambda)
{
  double x1old, x2old, dxold, fold, slope, F1, F2, tmplam, rhs1, rhs2, a, b, disc, lambd2, f2;
  int count; 

  x1old = x1;
  x2old = x2;
  dxold = dx;
  fold  = f;
 
  if (dxold > Smax)
    {
      dx1   = dx1 * Smax /dxold;
      dx2   = dx2 * Smax /dxold;
      dxold = Smax;
    }

  slope = g1*dx1 + g2*dx2;

  if (slope >= 0.)
    {
      printf ("Utility::NewtonBackTrack: Error - roundoff problem\n");
    }

  lambda  = 1.;

  for (int i = 0; i <= Maxiter; i++)
    {
      x1 = x1old + lambda * dx1;
      x2 = x2old + lambda * dx2;

      NewtonFunction (x1, x2, F1, F2);

      f = 0.5 * (F1*F1 + F2*F2);

      if (f <= fold + alpha * slope)
	{
	  dx = lambda * dxold;
	  
	  return;
	}
      else
	{
	  if (lambda == 1.)
	    tmplam = - slope /2. /(f - fold - slope);
	  else
	    {
	      rhs1 = f  - fold - lambda * slope;
	      rhs2 = f2 - fold - lambd2 * slope;

	      a = (           rhs1 /lambda/lambda -          rhs2 /lambd2/lambd2) /(lambda - lambd2);
	      b = (- lambd2 * rhs1 /lambda/lambda + lambda * rhs2 /lambd2/lambd2) /(lambda - lambd2);

	      if (a == 0.)
		tmplam = - slope /2./b;
	      else
		{
		  disc = b*b - 3. * a * slope;

		  if (disc < 0.)
		    tmplam = 0.5 * lambda;
		  else if (b <= 0.)
		    tmplam = (- b + sqrt (disc)) /3./a;
		  else
		    tmplam = - slope /(b + sqrt (disc));
		}

	      if (tmplam > 0.5 * lambda)
		tmplam = 0.5 * lambda;
	    }
	}

      lambd2 = lambda;
      f2     = f;
      lambda = Fmax (tmplam, 0.1*lambda);
     }

  lambda = Fmax (tmplam, 0.1*lambda);
  dx     = lambda * dxold;
	  
  return;
}

// #####################################################################
// Function to calculate Jacobian matrix for Newton-Raphson root finding
// #####################################################################
void Utility::NewtonJacobian (double x1, double x2, double& J11, double& J12, double& J21, double& J22)
{
  double F1m, F2m, F1p, F2p;
  
  NewtonFunction (x1 - dS, x2, F1m, F2m);
  NewtonFunction (x1 + dS, x2, F1p, F2p);

  J11 = (F1p - F1m) /2./dS;
  J21 = (F2p - F2m) /2./dS;

  NewtonFunction (x1, x2 - dS, F1m, F2m);
  NewtonFunction (x1, x2 + dS, F1p, F2p);

  J12 = (F1p - F1m) /2./dS;
  J22 = (F2p - F2m) /2./dS;
}

// ########################################
// Function to return maximum of two values
// ########################################
double Utility::Fmax (double f1, double f2)
{
  if (f1 > f2)
    return f1;
  else
    return f2;
}

// ########################################
// Function to return minimum of two values
// ########################################
double Utility::Fmin (double f1, double f2)
{
  if (f1 < f2)
    return f1;
  else
    return f2;
}

// ########################################
// Function to strip comments from a string
// ########################################
string Utility::StripComments (const string& input)
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
json Utility::ReadJSONFile (const string& filename)
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
	  JSONData = json::parse (StripComments (buffer.str ()));
        }
      catch (json::parse_error& e)
	{
	  cerr << "Utility::ReadJSONFile: Unable to parse JSON file: " << e.what() << endl;
	  exit (1);
        }
      JSONFile.close ();
    }
  else
    {
      cerr << "Utility::ReadJSONFile: Unable to open JSON file: " << filename << endl;
      exit (1);
    }

  return JSONData;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Utility::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Utility::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* Utility::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("Utility::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ################################################################
// Function to check that directory exists, and create it otherwise
// ################################################################
bool Utility::CreateDirectory (const char* path)
{
  struct stat st = {0};
  
  if (stat (path, &st) == -1)
    {
#ifdef _WIN32
      if (mkdir (path) != 0)
	{
	  printf ("Utility::CreateDirectory: Error creating directory: %s\n", path);
	  return false;
	}
#else
      if (mkdir (path, 0700) != 0)
	{
	  printf ("Utility::CreateDirectory: Error creating directory: %s\n", path);
	  return false;
	}
#endif
    }
  
  return true;
}

// #################################
// Function to call operating system
// #################################
void Utility::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("Utility::CallSystem: Operating system call error executing %s\n", command);
      exit (1);
    }
}

