// TearX.cpp

#include "TearX.h"

// ###########
// Constructor
// ###########
TearX::TearX ()
{
  // -------------------------------------------
  // Ensure that directory ../Outputs/Tear exits
  // -------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Tear"))
    {
      exit (1);
    }

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/TearX.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR   = JSONData["NTOR"]  .get<int>    ();
  q0     = JSONData["q0"]    .get<double> ();
  nu     = JSONData["nu"]    .get<double> ();
  alpha  = JSONData["alpha"] .get<double> ();

  eps    = JSONData["eps"]   .get<double> ();
  del    = JSONData["del"]   .get<double> ();
  EPS    = JSONData["EPS"]   .get<double> ();
  Psimax = JSONData["Psimax"].get<double> ();
  Nr     = JSONData["Nr"]    .get<int>    ();

  q0_sta = JSONData["q0_sta"].get<double> ();
  q0_end = JSONData["q0_end"].get<double> ();
  q0_num = JSONData["q0_num"].get<int>    ();

  nu_sta = JSONData["nu_sta"].get<double> ();
  nu_end = JSONData["nu_end"].get<double> ();
  nu_num = JSONData["nu_num"].get<int>    ();
  
  acc    = JSONData["acc"]   .get<double> ();
  h0     = JSONData["h0"]    .get<double> ();
  hmin   = JSONData["hmin"]  .get<double> ();
  hmax   = JSONData["hmax"]  .get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (q0 < 0.)
    {
      printf ("TearX:: Error - q0 cannot be negative\n");
      exit (1);
    }
  if (nu < 2.)
    {
      printf ("TearX:: Error - nu cannot be less than 2\n");
      exit (1);
    }
  if (alpha < 0.)
    {
      printf ("TearX:: Error - alpha cannot be negative\n");
      exit (1);
    }
 if (NTOR < 1)
    {
      printf ("TearX:: Error - NTOR must be positive\n");
      exit (1);
    }
  if (eps <= 0.)
    {
      printf ("TearX:: Error - eps must be positive\n");
      exit (1);
    }
  if (del <= 0.)
    {
      printf ("TearX:: Error - del must be positive\n");
      exit (1);
    }
  if (EPS <= 0.)
    {
      printf ("TearX:: Error - EPS must be positive\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("TearX:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (Psimax < 0. || Psimax > 1.)
    {
      printf ("TearX:: Error - Psimax must lie between 0 and 1\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("TearX:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("TearX:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("TearX:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("TearX:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (q0_num < 2)
    {
      printf ("TearX:: Error - q0_num must be greater than 1\n");
      exit (1);
    }
  if (nu_num < 2)
    {
      printf ("TearX:: Error - nu_num must be greater than 1\n");
      exit (1);
    }

  ntor = double (NTOR);
  qc   = nu * q0;

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  
  printf ("\nClass TEARX::\n\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ntor   = %-2d         q0     = %-10.3e nu     = %-10.3e qc  = %-10.3e alpha  = %-10.3e\n",
	  NTOR, q0, nu, qc, alpha);
  printf ("Nr     = %-5d      eps    = %-10.3e del    = %-10.3e EPS = %-10.3e Psimax = %-10.3e\n",
	  Nr, eps, del, EPS, Psimax);
  printf ("q0_sta = %-10.3e q0_end = %-10.3e q0_num = %-5d\n",
	  q0_sta, q0_end, q0_num);
  printf ("nu_sta = %-10.3e nu_end = %-10.3e nu_num = %-5d\n",
	  q0_sta, q0_end, q0_num);
  printf ("acc    = %-10.3e h0     = %-10.3e hmax   = %-10.3e\n",
	  acc, h0, hmax);
}

// ##########
// Destructor
// ##########
TearX::~TearX ()
{
}

// ##################
// Function to set q0
// ##################
void TearX::Setq0 (double _q0)
{
  q0 = _q0;
  qc = nu * q0;
}

// ##################
// Function to set nu
// ##################
void TearX::Setnu (double _nu)
{
  nu = _nu;
  qc = nu * q0;
}

// ###################
// Function to scan q0
// ###################
void TearX::Scanq0 ()
{
  FILE* file = OpenFilew ("../Outputs/TearX/Scanq0.txt");
  fclose (file);

  for (int i = 0; i < q0_num; i++)
    {
      double _q0 = q0_sta + (q0_end - q0_sta) * double (i) /double (q0_num - 1);

      Setq0 (_q0);

      Solve (1);
    }
}

// ###################
// Function to scan nu
// ###################
void TearX::Scannu ()
{
  FILE* file = OpenFilew ("../Outputs/TearX/Scannu.txt");
  fclose (file);

  for (int i = 0; i < nu_num; i++)
    {
      double _nu = q0_sta + (nu_end - nu_sta) * double (i) /double (nu_num - 1);

      Setnu (_nu);

      Solve (2);
    }
}

// #########################
// Function to solve problem
// #########################
void TearX::Solve (int flag)
{
  // ---------------
  // Allocate memory
  // ---------------
  rr  = new double[Nr];
  qq  = new double[Nr];
  PSI = new double[Nr];
  rho = new double[Nr];

  // ------------------
  // Set up radial grid
  // ------------------
  for (int i = 0; i < Nr; i++)
    {
      rr[i] = double (i) /double (Nr-1);
      qq[i] = Getq (rr[i]);
    }
  
  // -------------
  // Calculate PSI
  // -------------
  double r, h, t_err;
  int    rept;
  double y[1], err[1];

  rhs_chooser = 0;
  r           = eps;
  y[0]        = 0.;
  h           = h0;
  count       = 0;
  PSI[0]      = 0.;
  rho[0]      = 0.;

  for (int i = 1; i < Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (1, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (1, r, y, err, rr[i] - r);

      PSI[i] = y[0];
    }

  double PSIa = PSI[Nr-1];
  for (int i = 0; i < Nr; i++)
    {
      PSI[i] /= PSIa;
      rho[i] = sqrt (PSI[i]);
    }

  q_spline = gsl_spline_alloc (gsl_interp_cspline, Nr);
  P_spline = gsl_spline_alloc (gsl_interp_cspline, Nr);
  R_spline = gsl_spline_alloc (gsl_interp_cspline, Nr);

  q_acc = gsl_interp_accel_alloc ();
  P_acc = gsl_interp_accel_alloc ();
  R_acc = gsl_interp_accel_alloc ();

  gsl_spline_init (q_spline, rr,  qq,  Nr);
  gsl_spline_init (P_spline, PSI, rr,  Nr);
  gsl_spline_init (R_spline, rr,  rho, Nr);

  r95 = gsl_spline_eval (P_spline, 0.95, P_acc);
  R95 = gsl_spline_eval (R_spline, r95,  R_acc);
  q95 = Getq (r95);
  s95 = GetShear (r95);
  S95 = R95 * gsl_spline_eval_deriv (q_spline, r95, q_acc) /q95 /gsl_spline_eval_deriv (R_spline, r95, R_acc);

  rmax = gsl_spline_eval (P_spline, Psimax, P_acc);
  qmax = gsl_spline_eval (q_spline, rmax,   P_acc);

  if (flag == 0)
    printf ("\nr95  = %-10.3e R95  = %-10.3e q95  = %-10.3e s95 = %-10.3e S95 = %10.3e\n",
	    r95, R95, q95, s95, S95);
  if (flag == 0)
    printf ("rmax = %-10.3e qmax = %-10.3e\n\n",
	    rmax, qmax);

  // -------------------------------------
  // Determine number of rational surfaces
  // -------------------------------------
  nres = int (qmax /ntor) - int (q0 /ntor);
  if (nres <= 0)
    {
      printf ("No rational surfaces in plasma\n");
      exit (1);
    }

  // ---------------
  // Allocate memory
  // ---------------
  mres  = new int   [nres];
  rres  = new double[nres];
  Pres  = new double[nres];
  sres  = new double[nres];
  Dres  = new double[nres];
  Psi.resize (Nr, nres);
    
  ss    = new double[Nr];
  JJ    = new double[Nr];
  JJp   = new double[Nr];
  lvals = new double[Nr];
  
  // --------------------------------
  // Calculate equilibrium quantities
  // --------------------------------
  for (int i = 0; i < Nr; i++)
    {
      double q, s, J, Jp, lambda;
      
      GetEquilibrium (rr[i], q, s, J, Jp, lambda);

      ss[i]  = s;
      JJ[i]  = J;
      JJp[i] = Jp;
      
      if (i == 0)
	lvals[i] = 4.;
      else
	lvals[i] = rr[i]*lambda;
    }

  int mmax = int (qmax /ntor);
  int mmin = mmax - nres + 1;

  for (int isurf = 0; isurf < nres; isurf++)
    {
      // ----------------------
      // Find rational surfaces
      // ----------------------
      mpol = double (mmin + isurf);
      qs   = mpol /ntor;
      rs   = FindRationalSurface ();
      sr   = GetShear (rs);
 
      // ---------------------------------
      // Calculate tearing stability index
      // ---------------------------------
      Delta = GetDelta (isurf);
      if (flag == 0)
	printf ("mpol = %-2d ntor = %-2d rs = %-10.3e s = %-10.3e Delta = %11.4e\n",
		mmin + isurf, NTOR, rs, sr, Delta);

      // ----------
      // Store data
      // ----------
      mres[isurf] = mmin + isurf;
      rres[isurf] = rs;
      Pres[isurf] = pow (gsl_spline_eval (R_spline, rres[isurf],  R_acc), 2.);
      sres[isurf] = sr;
      Dres[isurf] = Delta;
    }

  // -----------------
  // Write netcdf file
  // -----------------
  if (flag == 0)
    WriteNetcdf ();

  // ---------------
  // Write scan data
  // ---------------
  if (flag == 1)
    {
      FILE* file = OpenFilea ("../Outputs/TearX/Scanq0.txt");
      
      printf ("q0 = %-10.3e qa/qc = %-10.3e q95 = %-10.3e nres = %-3d m = (%-3d, %-3d)\n",
	      q0, qq[Nr-1]/qc, q95, nres, mres[0], mres[nres-1]);
	      
      for (int i = 0; i < nres; i++)
	fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %3d %11.4e %11.4e %11.4e %11.4e\n",
		 q0, nu, qc, r95, q95, mres[i], rres[i], sres[i], Dres[i], Pres[i]);
      fclose (file);
    }
  if (flag == 2)
    {
      FILE* file = OpenFilea ("../Outputs/TearX/Scannu.txt");
      
      printf ("nu = %-10.3e qa/qc = %-10.3e q95 = %-10.3e nres = %-3d m = (%d, %d)\n",
	      nu, qq[Nr-1]/qc, q95, nres, mres[0], mres[nres-1]);
	      
      for (int i = 0; i < nres; i++)
	fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %3d %11.4e %11.4e %11.4e %11.4e\n",
		 q0, nu, qc, r95, q95, mres[i], rres[i], sres[i], Dres[i], Pres[i]);
      fclose (file);
    }
  
  // -------
  // Cleanup
  // -------
  delete[] rr;   delete[] ss;   delete[] JJ;   delete[] JJp; delete[] lvals;
  delete[] mres; delete[] rres; delete[] Dres; delete[] PSI; delete[] sres;
  delete[] rho;  delete[] Pres;

  gsl_spline_free (q_spline);    gsl_spline_free (P_spline);    gsl_spline_free (R_spline);
  gsl_interp_accel_free (q_acc); gsl_interp_accel_free (P_acc); gsl_interp_accel_free (R_acc);
}

// #####################################
// Function to write data to netcdf file
// #####################################
void TearX::WriteNetcdf ()
{
  double* Psi_y = new double[Nr*nres];

  int cnt = 0;
  for (int j = 0; j < Nr; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Psi_y[cnt] = Psi(j, jp);
	cnt++;
      }
  
  try
    {
      NcFile dataFile ("../Outputs/TearX/TearX.nc", NcFile::replace);
       
      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim i_d = dataFile.addDim ("Ni", 1);
      NcDim r_d = dataFile.addDim ("Nr", Nr);
      NcDim s_d = dataFile.addDim ("Ns", nres);

      vector<NcDim> psi_d;
      psi_d.push_back (r_d);
      psi_d.push_back (s_d);

      NcVar r_x   = dataFile.addVar ("r",      ncDouble, r_d);
      r_x.putVar (rr);
      NcVar q_x   = dataFile.addVar ("q",      ncDouble, r_d);
      q_x.putVar (qq);
      NcVar s_x   = dataFile.addVar ("s",      ncDouble, r_d);
      s_x.putVar (ss);
      NcVar J_x   = dataFile.addVar ("J",      ncDouble, r_d);
      J_x.putVar (JJ);
      NcVar Jp_x  = dataFile.addVar ("Jp",     ncDouble, r_d);
      Jp_x.putVar (JJp);
      NcVar l_x   = dataFile.addVar ("lambda", ncDouble, r_d);
      l_x.putVar (lvals);
      NcVar m_x   = dataFile.addVar ("mres",   ncInt, s_d);
      m_x.putVar (mres);
      NcVar rs_x  = dataFile.addVar ("rres",   ncDouble, s_d);
      rs_x.putVar (rres);
      NcVar Ps_x  = dataFile.addVar ("Pres",   ncDouble, s_d);
      Ps_x.putVar (Pres);
      NcVar ss_x  = dataFile.addVar ("sres",   ncDouble, s_d);
      ss_x.putVar (sres);
      NcVar D_x   = dataFile.addVar ("Dres",   ncDouble, s_d);
      D_x.putVar (Dres);
      NcVar psi_x = dataFile.addVar ("psi",    ncDouble, psi_d);
      psi_x.putVar (Psi_y);
    }
      catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/TearX/TearX.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] Psi_y;
}

// #########################################
// Function to return equilibrium quantities
// #########################################
void TearX::GetEquilibrium (double r, double& q, double& s, double& J, double& Jp, double& lambda)
{
  double r2  = r*r;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;

  q      = Getq (r);
  s      = r * Getqp (r) /q;
  J      = (2./q0) * pow (or2, nu1);
  Jp     = - (4./q0) * nu1 * r * pow (or2, nu2);
  if (r > 1.)
    {
      J  = 0.;
      Jp = 0.;
    }
  lambda = - q * Jp /s;
}


// #########################################
// Function to return value of safety-factor
// #########################################
double TearX::Getq (double r)
{
  double r2  = r*r;
  double r4  = r2*r2;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;
  double q;

  double f = (1. - pow (or2, nu)) /q0/nu;
  if (r > 1.)
    {
      f = 1. /q0/nu;
    }
  double lg = log (fabs (1. - r2) + EPS);
  
  if (r < 0.01)
    {
      q = q0 * (1. + 0.5*nu1*r2 + (0.25*nu1*nu1 - nu1*nu2/6.)*r4) + alpha*r2;
     }
  else
    {
      q = r2 /f - alpha * lg;
    }

  return q;
}

// ##############################################
// Function to return derivative of safety-factor
// ##############################################
double TearX::Getqp (double r)
{
  double r2  = r*r;
  double r3  = r*r2;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;
  double qp;

  double f  = (1. - pow (or2, nu)) /q0/nu;
  double fp = 2. * r * pow (or2, nu1) /q0;
  double lg = fabs (1. - r2) + EPS;

  if (r > 1.)
    {
      f  = 1. /q0/nu;
      fp = 0.;
    }
  
  if (r < 0.01)
    {
      qp = q0 * (nu1*r + 4.*(0.25*nu1*nu1 - nu1*nu2/6.)*r3) + 2.*alpha*r;
     }
  else
    {
      qp = 2.*r /f - r2 * fp /f/f + 2.*alpha*r /lg;
    }

  return qp;
}

// ##########################################
// Function to return value of magnetic shear
// ##########################################
double TearX::GetShear (double r)
{
  double q  = Getq (r);
  double qp = Getqp (r);

  return r * qp /q;
}
 
// ########################################
// Function to find rational surface radius
// ########################################
double TearX::FindRationalSurface ()
{
  double rs = -1.;
  
  int N = 10000;
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
	      printf ("TearX:FindRationalSurface - Error: rn is NaN\n");
	      exit (1);
	    }
	  else
	    rs = rn;

	  break;
	}
    }

  return rs;
}

// #############################################
// Function to calculate tearing stability index
// #############################################
double TearX::GetDelta (int isurf)
{
  // --------------------------------------------------------------------
  // Launch solution from magnetic axis and integrate to rational surface
  // --------------------------------------------------------------------
  double r, h, t_err;
  int    rept;
  double y[2], err[2];

  rhs_chooser = 1;
  r           = eps;
  y[0]        = pow (r, mpol);
  y[1]        = mpol * pow (r, mpol-1.);
  h           = h0;
  count       = 0;

  Psi(0, isurf) = y[0];

  int isave;
  for (int i = 1; i < Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);
      Psi(i, isurf)  = y[0];

      if (rr[i+1] > rs + del)
	{
	  isave = i;
	  break;
	}
    }

  do
    {
      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - del);
  CashKarp45Fixed (2, r, y, err, rs - del - r);

  // ---------------------------------------------------
  // Calculate coefficients of large and small solutions
  // ---------------------------------------------------
  double q, s, J, Jp, lambda;
  GetEquilibrium (rs, q, s, J, Jp, lambda);
  
  double a11    = 1. - lambda * del * (log (del) - 1.);
  double a12    =  - del;
  double a21    = lambda * log (del);
  double a22    = 1.;
  double jac    = a11 * a22 - a12 * a21;

  double Clm    = (  a22 * y[0] - a12 * y[1]) /jac;
  double Csm    = (- a21 * y[0] + a11 * y[1]) /jac;

  for (int i = 0; i <= isave; i++)
    Psi(i, isurf) /= Clm /mpol;
  
  // ----------------------------------------------------------------------
  // Launch solution from plasma boundary and integrate to rational surface
  // ----------------------------------------------------------------------
  r     = 1.;
  y[0]  = 1.;
  y[1]  = - mpol;
  h     = - h0;
  count = 0;

  Psi(Nr-1, isurf) = 1.;

  for (int i = Nr-2; i > 0; i--)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);
      Psi(i, isurf) = y[0];

      if (rr[i-1] < rs + del)
	{
	  isave = i;
	  break;
	}
    }
  
  do
    {
      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + del);
  CashKarp45Fixed (2, r, y, err, rs + del - r);

  a11 = 1. + lambda * del * (log (del) - 1.);
  a12 = del;
  a21 = lambda * log (del);
  a22 = 1.;
  jac = a11 * a22 - a12 * a21;

  double Clp = (  a22 * y[0] - a12 * y[1]) /jac;
  double Csp = (- a21 * y[0] + a11 * y[1]) /jac;

  for (int i = isave; i < Nr; i++)
    Psi(i, isurf) /= Clp /mpol;

  // ---------------------------------
  // Calculate tearing stability index
  // ---------------------------------
  double Delta_ = rs * (Csp /Clp - Csm /Clm);

  return Delta_;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void TearX::CashKarp45Rhs (double r, double* y, double* dydr)
{
  if (rhs_chooser == 0)
    {
      // y[0] = dPSI/dr
      
      double q = Getq (r);

      dydr[0] = r /q;
    }
  else if (rhs_chooser == 1)
    {
      // y[0] = psi
      // y[1] = dpsi/dr
  
      double q, s, J, Jp, lambda;
      
      GetEquilibrium (r, q, s, J, Jp, lambda);
      
      dydr[0] = y[1];
      dydr[1] = - y[1] /r + mpol*mpol * y[0] /r/r + Jp * y[0] /r /(1./q - 1./qs);
    }
  else
    {
      printf ("TearX: rhs_chooser error\n");
      exit (1);
    }
 }

