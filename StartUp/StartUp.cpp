// StartUp.cpp

#include "StartUp.h"

#define NINPUT 1

// ###########
// Constructor
// ###########
StartUp::StartUp ()
{
  // --------------------------------------
  // Ensure that directory ../Outputs exits
  // --------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/StartUp.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  alpha = JSONData["alpha"].get<double> ();
  eps   = JSONData["eps"]  .get<double> ();
  Nx    = JSONData["Nx"]   .get<int>    ();
  
  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  Nint    = JSONData["Nint"]   .get<int>    ();
  Eta     = JSONData["Eta"]    .get<double> ();
  Maxiter = JSONData["Maxiter"].get<int>    ();
  lmin    = JSONData["lmin"]   .get<double> ();
  lmax    = JSONData["lmax"]   .get<double> ();
  xmax    = JSONData["xmax"]   .get<double> ();
  amax    = JSONData["amax"]   .get<double> ();
  Na      = JSONData["Na"]     .get<int>    ();
  
  // ---------------------
  // Print welcome message
  // ---------------------
  printf ("\nClass StartUp::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");

   // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("Calculation Parameters:\n");
  printf ("alpha = %11.4e eps  = %11.4e Nx      =  %4d\n",
	  alpha, eps, Nx);
  printf ("acc   = %11.4e h0   = %11.4e hmin    = %11.4e hmax = %11.4e\n",
	  acc, h0, hmin, hmax);
  printf ("Nint  = %3d         Eta  = %11.4e Maxiter = %3d         lmin = %11.4e lmax = %11.4e\n",
	  Nint, Eta, Maxiter, lmin, lmax);
  printf ("xmax  = %11.4e amax = %11.4e Na      =  %4d\n",
	  xmax, amax, Na);
}

// ##########
// Destructor
// ##########
StartUp::~StartUp ()
{
  delete[] xx;  delete[] T;   delete[] Bt;  delete[] q;
  delete[] aa;  delete[] qqa; delete[] lli; delete[] TTr;
  delete[] ttr; delete[] EEr; delete[] ll;
}

// #########################
// Function to solve problem
// #########################
void StartUp::Solve ()
{
  // ----------------
  // Calculate lambda
  // ----------------
  double _lambda = RootFind (lmin, lmax);

  // ------------------
  // Calculate profiles
  // ------------------
  xx = new double[Nx];
  T  = new double[Nx];
  Bt = new double[Nx];
  q  = new double[Nx];

  for (int i = 0; i < Nx; i++)
    xx[i] = double (i) /double (Nx-1);

  T[0]  = 1.;
  Bt[0] = 0.;

  double  x, h, t_err;
  int     rept; rhs_chooser = 1;
  double* y   = new double[5];
  double* err = new double[5];
 
  count = 0;
  h     = h0;
  x     = eps;
  y[0]  = 1. - lambda * eps*eps /2.;
  y[1]  = - lambda * eps;
  y[2]  = eps /2.;
  y[3]  = 0.;
  y[4]  = 0.;

  T [0] = 1.;
  Bt[0] = 0.; 
  for (int i = 1; i < Nx; i++)
    {
      do
	{
	  CashKarp45Adaptive (5, x, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (x < xx[i]);
      CashKarp45Fixed (5, x, y, err, xx[i] - x);

      T [i] = y[0];
      Bt[i] = y[2];
    }

  Bt1   = y[2];
  qa    = 1. /2. /Bt1;
  li    = 2. * y[4] /Bt1/Bt1;
  Tramp = 1. /pow (lambda, 0.4) /pow (Bt1, 0.8);
  Eramp = pow (lambda*lambda*lambda * Bt1, 0.2);
  tramp = y[3] /Bt1 /2. /pow (lambda*lambda*lambda * Bt1, 0.2);

  q[0] = 1.;
  for (int i = 0; i < Nx; i++)
    {
      q[i] = xx[i] /2./ Bt[i];
    }
  
  delete[] y; delete[] err;
  
  printf ("\nalpha = %11.4e lambda = %11.4e Bt1   = %11.4e qa = %11.4e  li   = %11.4e\n",
	  alpha, _lambda, Bt1, qa, li);
  printf ("Tramp = %11.4e Eramp  = %11.4e tramp = %11.4e\n\n",
	  Tramp, Eramp, tramp);

  // ------------------
  // Perform alpha scan
  // ------------------
  aa  = new double[Na];
  ll  = new double[Na];
  qqa = new double[Na];
  lli = new double[Na];
  TTr = new double[Na];
  ttr = new double[Na];
  EEr = new double[Na];

  for (int i = 0; i < Na; i++)
    {
      alpha = amax * double (i) /double (Na - 1);

      SolveScan ();
      if (i%100 == 0)
	printf ("alpha = %11.4e lambda = %11.4e qa = %11.4e li = %11.4e\n",
		alpha, lambda, qa, li);

      aa [i] = alpha;
      ll [i] = lambda;
      qqa[i] = qa;
      lli[i] = li;
      TTr[i] = Tramp;
      ttr[i] = tramp;
      EEr[i] = Eramp;
    }
  printf ("\n");
   
  // -----------------
  // Write netcdf file
  // -----------------
  Write_netcdf ();
}

// ########################################
// Function to solve problem for alpha scan
// ########################################
void StartUp::SolveScan ()
{
  double _lambda = RootFind (lmin, lmax);

  double  x, h, t_err;
  int     rept; rhs_chooser = 1;
  double* y   = new double[5];
  double* err = new double[5];
 
  count = 0;
  h     = h0;
  x     = eps;
  y[0]  = 1. - lambda * eps*eps /2.;
  y[1]  = - lambda * eps;
  y[2]  = eps /2.;
  y[3]  = 0.;
  y[4]  = 0.;

  do
    {
      CashKarp45Adaptive (5, x, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (x < 1.);
  CashKarp45Fixed (5, x, y, err, 1. - x);

  Bt1   = y[2];
  qa    = 1. /2. /Bt1;
  li    = 2. * y[4] /Bt1/Bt1;
  Tramp = 1. /pow (lambda, 0.4) /pow (Bt1, 0.8);
  Eramp = pow (lambda*lambda*lambda * Bt1, 0.2);
  tramp = y[3] /Bt1 /2. /pow (lambda*lambda*lambda * Bt1, 0.2);

  delete[] y; delete[] err;
} 

// ###################################
// Right hand sides of layer equations
// ###################################
void StartUp::CashKarp45Rhs (double x, double* y, double* dydx)
{
  if (rhs_chooser == 0)
    {
      // y[0] = T
      // y[1] = dT/dx

      double fa = (1. + alpha) /(pow (2., alpha + 1.) - 1.);
      double f1 = 2. * alpha * x /(1. + x*x);
      double f2 = fa * pow (1. + x*x, alpha);

      double f3;
      if (y[0] > 0.)
	f3 =   pow (  y[0], 1.5);
      else
	f3 = - pow (- y[0], 1.5);

      dydx[0] =   y[1];
      dydx[1] = - y[1] /x - f1 * y[1] - lambda * f3 /f2;
     }
  else
    {
      // y[0] = T
      // y[1] = dT/dx
      // y[2] = Bt
      // y[3] = int Bt
      // y[4] = int x Bt^2

      double fa = (1. + alpha) /(pow (2., alpha + 1.) - 1.);
      double f1 = 2. * alpha * x /(1. + x*x);
      double f2 = fa * pow (1. + x*x, alpha);

      double f3;
      if (y[0] > 0.)
	f3 =   pow (  y[0], 1.5);
      else
	f3 = - pow (- y[0], 1.5);

      dydx[0] =   y[1];
      dydx[1] = - y[1] /x - f1 * y[1] - lambda * f3 /f2;
      dydx[2] = - y[2] /x + f3;
      dydx[3] =   y[2];
      dydx[4] = x * y[2]*y[2];
    }
}
  
// ################################
// Target function for root finding
// ################################
double StartUp::RootFindF (double _lambda)
{
  lambda = _lambda;
  
  double  x, h, t_err;
  int     rept; rhs_chooser = 0;
  double* y = new double[2];

  count = 0;
  h     = h0;
  x     = eps;
  y[0]  = 1. - lambda * eps*eps /2.;
  y[1]  = - lambda * eps;

  double xold, yold;
  do
    {
      xold = x;
      yold = y[0];
      
      CashKarp45Adaptive (2, x, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (y[0] > 0. && x < xmax);

  double xzero;
  if (x > xmax)
    xzero = xmax;
  else
    xzero = xold + yold * (x - xold) /(yold - y[0]);

  //printf ("lambda = %11.4e xzero = %11.4e\n", lambda, xzero);

  delete[] y;
  
  return xzero - 1.;
}

// #############################
// Function to write netcdf file
// #############################
void StartUp::Write_netcdf ()
{
  printf ("Writing data to netcdf file Outputs/StartUp.nc:\n");

   double Input[NINPUT];
   
   Input[0] = alpha;
   
   try
     {
       NcFile dataFile ("../Outputs/StartUp.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim i_d = dataFile.addDim ("Ni", NINPUT);
       NcDim x_d = dataFile.addDim ("Nx", Nx);
       NcDim a_d = dataFile.addDim ("Na", Na);
 
       NcVar i_x = dataFile.addVar ("InputParameters", ncDouble, i_d);
       i_x.putVar (Input);
       NcVar x_x = dataFile.addVar ("x",               ncDouble, x_d);
       x_x.putVar (xx);
       NcVar T_x = dataFile.addVar ("T_e",             ncDouble, x_d);
       T_x.putVar (T);
       NcVar q_x = dataFile.addVar ("q",               ncDouble, x_d);
       q_x.putVar (q);
       NcVar a_x = dataFile.addVar ("alpha",           ncDouble, a_d);
       a_x.putVar (aa);
       NcVar z_x = dataFile.addVar ("lambda",          ncDouble, a_d);
       z_x.putVar (ll);
       NcVar p_x = dataFile.addVar ("q_a",             ncDouble, a_d);
       p_x.putVar (qqa);
       NcVar l_x = dataFile.addVar ("l_i",             ncDouble, a_d);
       l_x.putVar (lli);
       NcVar Q_x = dataFile.addVar ("T_ramp",          ncDouble, a_d);
       Q_x.putVar (TTr);
       NcVar t_x = dataFile.addVar ("t_ramp",          ncDouble, a_d);
       t_x.putVar (ttr);
       NcVar E_x = dataFile.addVar ("E_ramp",          ncDouble, a_d);
       E_x.putVar (EEr);
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/ECE/ECE.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }
}
