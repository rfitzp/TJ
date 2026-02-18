// Kink.cpp

#include "Kink.h"

// ###########
// Constructor
// ###########
Kink::Kink ()
{
  // -------------------------------------------
  // Ensure that directory ../Outputs/Kink exits
  // -------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Kink"))
    {
      exit (1);
    }

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/Kink.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR = JSONData["NTOR"].get<int>    ();
  MPOL = JSONData["MPOL"].get<int>    ();
  qa   = JSONData["qa"]  .get<double> ();
  nu   = JSONData["nu"]  .get<double> ();

  Nqa  = JSONData["Nqa"].get<int> ();
  Nnu  = JSONData["Nnu"].get<int> ();

  eps  = JSONData["eps"].get<double> ();
  del  = JSONData["del"].get<double> ();
  Nr   = JSONData["Nr"] .get<int>    ();

  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  Eta     = JSONData["Eta"]    .get<double> ();
  Maxiter = JSONData["Maxiter"].get<int>    ();
  Nint    = JSONData["Nint"]   .get<int>    ();

  // ------------
  // Sanity check
  // ------------
  if (qa < 0.)
    {
      printf ("Kink:: Error - qa cannot be negative\n");
      exit (1);
    }
  if (nu < 1.)
    {
      printf ("Kink:: Error - nu cannot be less than 1\n");
      exit (1);
    }
 if (NTOR < 1)
    {
      printf ("Kink:: Error - NTOR must be positive\n");
      exit (1);
    }
 if (MPOL < 1)
    {
      printf ("Kink:: Error - MPOL must be positive\n");
      exit (1);
    }
  if (eps <= 0.)
    {
      printf ("Kink:: Error - eps must be positive\n");
      exit (1);
    }
  if (del <= 0.)
    {
      printf ("Kink:: Error - del must be positive\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("Kink:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("Kink:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("Kink:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("Kink:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("Kink:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("Kink:: Error - hmax must exceed hmin\n");
      exit (1);
    }

  ntor = double (NTOR);
  mpol = double (MPOL);

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\nClass KINK::\n\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ntor = %-2d   mpol = %-2d         qa = %9.3e nu  = %9.3e q0 = %9.3e\n",
	  NTOR, MPOL, qa, nu, qa /(1. + nu));
  printf ("Nr   = %4d eps  = %9.3e del = %9.3e acc = %9.3e h0 = %9.3e hmax = %9.3e\n",
	  Nr, eps, del, acc, h0, hmax);
}

// ##########
// Destructor
// ##########
Kink::~Kink ()
{
}

// #########################
// Function to solve problem
// #########################
void Kink::Solve ()
{
  // ---------------
  // Allocate memory
  // ---------------
  rr   = new double[Nr];
  qq   = new double[Nr];
  ss   = new double[Nr];
  ppsi = new double[Nr];

  qqa = new double[Nqa];
  ll  = new double[Nqa];
  gg  = new double[Nqa];
  bb  = new double[Nqa];
  ww  = new double[Nqa];

  nn  = new double[Nnu];
  qxa = new double[Nnu];

  // ------------------
  // Set up radial grid
  // ------------------
  for (int i = 0; i < Nr; i++)
    {
      rr[i] = double (i) /double (Nr-1);
      
      double q, sigma, Drive;

      GetEquilibrium (rr[i], q, sigma, Drive);

      qq[i]   = q;
      ss[i]   = sigma;
      ppsi[i] = 0.;
    }

  // ---------------------
  // Find rational surface
  // ---------------------
  qs = mpol /ntor;
  rs = FindRationalSurface ();
  
  if (rs > 0. && rs < 1.)
    {
      printf ("\nFound rational surface: mpol = %-2d ntor = %-2d rs = %10.4e residual = %10.4e\n",
	      MPOL, NTOR, rs, fabs (Getq(rs) - qs));
    }
  else
    {
      printf ("\nmpol = %-2d ntor = %-2d: No rational surface in plasma\n",
	      MPOL, NTOR);
    }

  // -----------------------------------
  // Calculate lambda, gamma, and b_crit
  // -----------------------------------
  double lambda = GetLambda ();
  double gamma  = - (lambda + fabs (mpol));
  double bcrit  = 0.;

  if (gamma > 0.)
    bcrit = pow ((lambda - fabs (mpol)) / (lambda + fabs (mpol)), 1./2. /fabs(mpol));

  printf ("\nlambda = %11.4e  gamma = %11.4e  b_crit = %11.4e\n\n", lambda, gamma, bcrit);

  // -------
  // Scan qa
  // -------
  //Scanqa ();

  // ---------------------------------------------
  // Find critical qa value for marginal stability
  // ---------------------------------------------
  FindCritical ();
  
  // -----------------
  // Write netcdf file
  // -----------------
  WriteNetcdf ();
  
  // -------
  // Cleanup
  // -------
  delete[] rr;  delete[] qq; delete[] ppsi;
  delete[] qqa; delete[] ll; delete[] gg;   delete[] bb; delete[] ww;
  delete[] qxa; delete[] nn;
}

// ###################################################
// Function to find critical qa for marginal stability
// ###################################################
void Kink::FindCritical ()
{
  // ----------------------------------------------
  // Scan nu to find critical value of qa versus nu
  // ----------------------------------------------
  for (int i = 0; i < Nnu; i++)
    {
      double nu_min = 1.001;
      double nu_max = 2.0;

      nu = nu_min + (nu_max - nu_min) * double (i) /double (Nnu - 1);
      
      double x;
      double x1 = 0.5    * (mpol /ntor);
      double x2 = 0.9999 * (mpol /ntor);
      double F1 = RootFindF (x1);
      double F2 = RootFindF (x2);

      //x = RootFind (x1, x2);
      
      //nn[i]  = nu;
      //qxa[i] = x;
      
      if (F1*F2 < 0.)
	{
	  Ridder (x1, x2, F1, F2, x);

	  nn[i]  = nu;
	  qxa[i] = x;
	}
      else
	{
	  nn[i]  = 1.;
	  qxa[i] = mpol/ntor;
	}
    }
}

// ###################
// Function to scan qa
// ###################
void Kink::Scanqa ()
{
  double qamin = 0.01    * (mpol /ntor);
  double qamax = 0.99999 * (mpol /ntor);
  
  for (int i = 0; i < Nqa; i++)
    {
      // -------
      // Scan qa
      // -------
      qa = qamin + (qamax - qamin) * double (i) /double (Nqa - 1);
      
      // ---------------------
      // Find rational surface
      // ---------------------
      qs = mpol /ntor;
      rs = FindRationalSurface ();
      
      // ---------------------------------------
      // Calculate lambda, gamma, b_crit, and dW
      // ---------------------------------------
      double lambda = GetLambda ();
      double gamma  = - (lambda + fabs (mpol));
      double bcrit  = 0.90;
      double dW     = (1./qa - 1./qs) * (1./qa - 1./qs) * (lambda + fabs (mpol));
      
      if (gamma > 0.)
	bcrit = pow ((lambda - fabs (mpol)) / (lambda + fabs (mpol)), 1./2. /fabs(mpol));

      qqa[i] = qa;
      ll[i]  = lambda;
      gg[i]  = gamma;
      bb[i]  = bcrit;
      if (dW < 0.)
	ww[i] = - dW;
      else
	ww[i] = - 1.e-1;
    }
}

// #####################################
// Function to write data to netcdf file
// #####################################
void Kink::WriteNetcdf ()
{
  double para[4];
  para[0] = rs;
  para[1] = mpol;
  para[2] = ntor;
  para[3] = nu;
  
  try
    {
      NcFile dataFile ("../Outputs/Kink/Kink.nc", NcFile::replace);
       
      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim p_d = dataFile.addDim ("Np",  4);
      NcDim r_d = dataFile.addDim ("Nr",  Nr);
      NcDim q_d = dataFile.addDim ("Nqa", Nqa);
      NcDim n_d = dataFile.addDim ("Nnu", Nnu);
      
      NcVar p_x = dataFile.addVar ("para", ncDouble, p_d);
      p_x.putVar (para);

      NcVar r_x   = dataFile.addVar ("r",     ncDouble, r_d);
      r_x.putVar (rr);
      NcVar q_x   = dataFile.addVar ("q",     ncDouble, r_d);
      q_x.putVar (qq);
      NcVar s_x   = dataFile.addVar ("sigma", ncDouble, r_d);
      s_x.putVar (ss);
      NcVar psi_x = dataFile.addVar ("psi",   ncDouble, r_d);
      psi_x.putVar (ppsi);

      NcVar qa_x = dataFile.addVar ("q_a",    ncDouble, q_d);
      qa_x.putVar (qqa);
      NcVar ll_x = dataFile.addVar ("lambda", ncDouble, q_d);
      ll_x.putVar (ll);
      NcVar gg_x = dataFile.addVar ("gamma",  ncDouble, q_d);
      gg_x.putVar (gg);
      NcVar bb_x = dataFile.addVar ("b_crit", ncDouble, q_d);
      bb_x.putVar (bb);
      NcVar ww_x = dataFile.addVar ("dW_nw",  ncDouble, q_d);
      ww_x.putVar (ww);

      NcVar qx_x = dataFile.addVar ("q_xa", ncDouble, n_d);
      qx_x.putVar (qxa);
      NcVar nn_x = dataFile.addVar ("nu",   ncDouble, n_d);
      nn_x.putVar (nn);
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Kink/Kink.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }
}

// #########################################
// Function to return equilibrium quantities
// #########################################
void Kink::GetEquilibrium (double r, double& q, double& sigma, double& Drive)
{
  double r2  = r*r;
  double r4  = r2*r2;
  double or2 = 1. - r2;
  double np1 = nu + 1.;
  double nm1 = nu - 1.;
  double q0  = qa /np1;

  double f = (1. - pow (or2, np1)) /qa;

  if (r > 1.)
    {
      q = qa * r*r;
    }
  else if (r < 0.01)
    {
      q = q0 * (1. + 0.5*nu*r2 + (0.25*nu*nu - nu*nm1/6.)*r4);
    }
  else
    {
      q = r2 /f;
    }

  if (r > 1.)
    {
      sigma = 0.;
      Drive = 0.;
    }
  else
    {
      sigma = (2./q0) * pow (1. - r2, nu);
      Drive = 4. * nu*np1 * pow (1. - r2, nm1) /qa;
    }
}

// #########################################
// Function to return value of safety-factor
// #########################################
double Kink::Getq (double r)
{
  double r2  = r*r;
  double r4  = r2*r2;
  double or2 = 1. - r2;
  double np1 = nu + 1.;
  double nm1 = nu - 1.;
  double q0  = qa /np1;
  double q;

  double f = (1. - pow (or2, np1)) /qa;
  if (r > 1.)
    {
      q = qa *r*r;
    }
  else if (r < 0.01)
    {
      q = q0 * (1. + 0.5*nu*r2 + (0.25*nu*nu - nu*nm1/6.)*r4);
     }
  else
    {
      q = r2 /f;
    }

  return q;
}

// ##############################################
// Function to return derivative of safety-factor
// ##############################################
double Kink::Getqp (double r)
{
  double r2  = r*r;
  double r3  = r*r2;
  double or2 = 1. - r2;
  double np1 = nu + 1.;
  double nm1 = nu - 1.;
  double q0  = qa /np1;
  double qp;

  double f  = (1. - pow (or2, np1))  /qa;
  double fp = 2. * r * pow (or2, nu) /q0;
  
  if (r < 0.01)
    {
      qp = q0 * (nu*r + 4.*(0.25*nu*nu - nu*nm1/6.)*r3);
     }
  else
    {
      qp = 2.*r /f - r2 * fp /f/f;
    }

  return qp;
}
 
// ########################################
// Function to find rational surface radius
// ########################################
double Kink::FindRationalSurface ()
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
	      printf ("Kink:FindRationalSurface - Error: rn is NaN\n");
	      exit (1);
	    }
	  else
	    rs = rn;

	  break;
	}
    }

  return rs;
}

// ##########################################
// Function to calculate kink stability index
// ##########################################
double Kink::GetLambda ()
{
  double psia, psipa;
  
  if (rs > 0. && rs < 1.)
    {
      // ----------------------------------------------------------------------
      // Launch solution from rational surface and integrate to plasma boundary
      // ----------------------------------------------------------------------
      double r, h, t_err;
      int    rept;
      double y[2], err[2];

      r     = rs + del;
      y[0]  = 0.;
      y[1]  = 1.;
      h     = h0;
      count = 0;

      for (int i = 1; i < Nr; i++)
	{
	  if (rr[i] > rs + del)
	    {
	      do
		{
		  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
		}
	      while (r < rr[i]);
	      CashKarp45Fixed (2, r, y, err, rr[i] - r);
	      ppsi[i] = y[0];
	    }
	}

      psia  = y[0];
      psipa = y[1];
    }
  else
    {
      // -------------------------------------------------------------------
      // Launch solution from magnetic axis and integrate to plasma boundary
      // -------------------------------------------------------------------
      double r, h, t_err;
      int    rept;
      double y[2], err[2];
      
      r     = eps;
      y[0]  = pow (r, mpol);
      y[1]  = mpol * pow (r, mpol-1.);
      h     = h0;
      count = 0;
      
      ppsi[0] = y[0];
      
      for (int i = 1; i < Nr; i++)
	{
	  do
	    {
	      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < rr[i]);
	  CashKarp45Fixed (2, r, y, err, rr[i] - r);
	  ppsi[i]  = y[0];
	}

      psia  = y[0];
      psipa = y[1];
    }

  // ------------------
  // Normalize solution
  // ------------------
  for (int i = 0; i < Nr; i++)
    ppsi [i] /= psia;

  // ----------------
  // Calculate lambda
  // ----------------
  return psipa /psia;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Kink::CashKarp45Rhs (double r, double* y, double* dydr)
{
  // y[0] = psi
  // y[1] = dpsi/dr
  
  double q, sigma, Drive;

  GetEquilibrium (r, q, sigma, Drive);

  dydr[0] = y[1];
  dydr[1] = - y[1] /r + mpol*mpol * y[0] /r/r - Drive * y[0] /(1./q - 1./qs); 
 }

// ################################
// Target function for zero finding
// ################################
double Kink::RootFindF (double x)
{
  qa = x;

  double lambda = GetLambda ();

  return lambda + fabs (mpol);
}
