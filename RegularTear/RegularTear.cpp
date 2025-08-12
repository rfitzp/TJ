// RegularTear.cpp

#include "RegularTear.h"

// ###########
// Constructor
// ###########
RegularTear::RegularTear ()
{
  // -------------------------------------------
  // Ensure that directory ../Outputs/RegularTear exits
  // -------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/RegularTear"))
    {
      exit (1);
    }

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/RegularTear.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR  = JSONData["NTOR"] .get<int>    ();
  q0    = JSONData["q0"]   .get<double> ();
  qa    = JSONData["qa"]   .get<double> ();
  eps   = JSONData["eps"]  .get<double> ();
  del   = JSONData["del"]  .get<double> ();
  Nr    = JSONData["Nr"]   .get<int>    ();
  acc   = JSONData["acc"]  .get<double> ();
  h0    = JSONData["h0"]   .get<double> ();
  hmin  = JSONData["hmin"] .get<double> ();
  hmax  = JSONData["hmax"] .get<double> ();
  Fixed = JSONData["Fixed"].get<int>    ();
  sigma = JSONData["sigma"].get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (q0 < 0.)
    {
      printf ("RegularTear:: Error - q0 cannot be negative\n");
      exit (1);
    }
  if (qa < 2.*q0)
    {
      printf ("RegularTear:: Error - qa cannot be less than 2*q0\n");
      exit (1);
    }
 if (NTOR < 1)
    {
      printf ("RegularTear:: Error - NTOR must be positive\n");
      exit (1);
    }
  if (eps <= 0.)
    {
      printf ("RegularTear:: Error - eps must be positive\n");
      exit (1);
    }
  if (del <= 0.)
    {
      printf ("RegularTear:: Error - del must be positive\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("RegularTear:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("RegularTear:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("RegularTear:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("RegularTear:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("RegularTear:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("RegularTear:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (sigma <= 0. || sigma > 1.)
    {
      printf ("RegularTear:: Error - sigma must lie in range 0 to 1\n");
      exit (1);
    }

  ntor = double (NTOR);
  nu   = qa /q0;

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  
  printf ("\nClass RegularTear::\n\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ntor = %-2d    q0  = %10.3e qa  = %10.3e sigma = %10.3e Fixed = %1d\n",
	  NTOR, q0, qa, sigma, Fixed);
  printf ("Nr   = %4d eps  = %10.3e del = %10.3e acc   = %10.3e h0    = %10.3e hmax = %10.3e\n",
	  Nr, eps, del, acc, h0, hmax);
}

// ##########
// Destructor
// ##########
RegularTear::~RegularTear ()
{
}

// #########################
// Function to solve problem
// #########################
void RegularTear::Solve ()
{
  // -------------------------------------
  // Determine number of rational surfaces
  // -------------------------------------
  nres = int (qa /ntor) - int (q0 /ntor);
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
  sres  = new double[nres];
  Dres  = new double[nres];
  Dresr = new double[nres];

  Psi  .resize (Nr, nres);
  Psip .resize (Nr, nres);
  Psipp.resize (Nr, nres);
  Cur  .resize (Nr, nres);

  Psir  .resize (Nr, nres);
  Psipr .resize (Nr, nres);
  Psippr.resize (Nr, nres);
  Curr  .resize (Nr, nres);

  Norm.resize (nres);
  Norp.resize (nres);
    
  rr    = new double[Nr];
  qq    = new double[Nr];
  ss    = new double[Nr];
  JJ    = new double[Nr];
  JJp   = new double[Nr];
  lvals = new double[Nr];
  
  // ------------------
  // Set up radial grid
  // ------------------
  for (int i = 0; i < Nr; i++)
    {
      rr[i] = double (i) /double (Nr-1);

      double q, s, J, Jp, lambda;
      
      GetEquilibrium (rr[i], q, s, J, Jp, lambda);

      qq[i]  = q;
      ss[i]  = s;
      JJ[i]  = J;
      JJp[i] = Jp;
      
      if (i == 0)
	lvals[i] = 4.;
      else
	lvals[i] = rr[i]*lambda;
    }

  int mmax = int (qa /ntor);
  int mmin = mmax - nres + 1;

  for (int isurf = 0; isurf < nres; isurf++)
    {
      // ----------------------
      // Find rational surfaces
      // ----------------------
      mpol = double (mmin + isurf);
      qs   = mpol /ntor;
      rs   = FindRationalSurface ();

      double q, J, Jp, lambda;
      
      GetEquilibrium (rs, q, sh, J, Jp, lambda);
  
      if (rs > 0. && rs < 1.)
	{
	  printf ("\nFound rational surface: mpol = %-2d ntor = %-2d rs = %11.4e ss = %11.4e residual = %11.4e\n",
		  mmin + isurf, NTOR, rs, sres, fabs (Getq(rs) - qs));
	}
      else
	{
	  printf ("\nmpol = %-2d ntor = %-2d: No rational surface in plasma\n",
		  mmin + isurf, NTOR);
	  exit (0);
	}

      // ---------------------------------
      // Calculate tearing stability index
      // ---------------------------------
      Delta = GetDelta (isurf);
      printf ("\nDelta = %11.4e\n", Delta);

      // --------------------------------------------
      // Calculate regularized tearing stablity index
      // --------------------------------------------
      Deltr = GetRegularDelta (isurf);
      printf ("\nDeltr = %11.4e\n", Deltr);

      // ----------
      // Store data
      // ----------
      mres[isurf] = mmin + isurf;
      rres[isurf] = rs;
      Dres[isurf] = Delta;
      sres[isurf] = sh;
    }
 
  // -----------------
  // Write netcdf file
  // -----------------
  WriteNetcdf ();
  
  // -------
  // Cleanup
  // -------
  delete[] rr;   delete[] ss;   delete[] JJ;   delete[] JJp;   delete[] lvals;
  delete[] mres; delete[] rres; delete[] Dres; delete[] Dresr; delete[] sres;
}

// #####################################
// Function to write data to netcdf file
// #####################################
void RegularTear::WriteNetcdf ()
{
  double  sig_y[1];
  double* Psi_y = new double[Nr*nres];
  double* Psp_y = new double[Nr*nres];
  double* Ppp_y = new double[Nr*nres];
  double* Cur_y = new double[Nr*nres];
  double* Psi_z = new double[Nr*nres];
  double* Psp_z = new double[Nr*nres];
  double* Ppp_z = new double[Nr*nres];
  double* Cur_z = new double[Nr*nres];

  sig_y[0] = sigma;
  
  int cnt = 0;
  for (int j = 0; j < Nr; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Psi_y[cnt] = Psi  (j, jp);
	Psp_y[cnt] = Psip (j, jp);
	Ppp_y[cnt] = Psipp(j, jp);
	Cur_y[cnt] = Cur  (j, jp);

	Psi_z[cnt] = Psir  (j, jp);
	Psp_z[cnt] = Psipr (j, jp);
	Ppp_z[cnt] = Psippr(j, jp);
	Cur_z[cnt] = Curr  (j, jp);

	cnt++;
      }
  
  try
    {
      NcFile dataFile ("../Outputs/RegularTear/RegularTear.nc", NcFile::replace);
       
      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim i_d = dataFile.addDim ("Ni", 1);
      NcDim r_d = dataFile.addDim ("Nr", Nr);
      NcDim s_d = dataFile.addDim ("Ns", nres);

      vector<NcDim> psi_d;
      psi_d.push_back (r_d);
      psi_d.push_back (s_d);

      NcVar sig_x  = dataFile.addVar ("sigma",   ncDouble, i_d);
      sig_x.putVar (sig_y);
      NcVar r_x    = dataFile.addVar ("r",       ncDouble, r_d);
      r_x.putVar (rr);
      NcVar q_x    = dataFile.addVar ("q",       ncDouble, r_d);
      q_x.putVar (qq);
      NcVar s_x    = dataFile.addVar ("s",       ncDouble, r_d);
      s_x.putVar (ss);
      NcVar J_x    = dataFile.addVar ("J",       ncDouble, r_d);
      J_x.putVar (JJ);
      NcVar Jp_x   = dataFile.addVar ("Jp",      ncDouble, r_d);
      Jp_x.putVar (JJp);
      NcVar l_x    = dataFile.addVar ("lambda",  ncDouble, r_d);
      l_x.putVar (lvals);
      NcVar m_x    = dataFile.addVar ("mres",    ncInt,    s_d);
      m_x.putVar (mres);
      NcVar rs_x   = dataFile.addVar ("rres",    ncDouble, s_d);
      rs_x.putVar (rres);
      NcVar ss_x   = dataFile.addVar ("sres",    ncDouble, s_d);
      ss_x.putVar (sres);
      NcVar D_x    = dataFile.addVar ("Dres",    ncDouble, s_d);
      D_x.putVar (Dres);
      NcVar psi_x  = dataFile.addVar ("psi",     ncDouble, psi_d);
      psi_x.putVar (Psi_y);
      NcVar psp_x  = dataFile.addVar ("psi_p",   ncDouble, psi_d);
      psp_x.putVar (Psp_y);
      NcVar ppp_x  = dataFile.addVar ("psi_pp",  ncDouble, psi_d);
      ppp_x.putVar (Ppp_y);
      NcVar cur_x  = dataFile.addVar ("j",       ncDouble, psi_d);
      cur_x.putVar (Cur_y);
      NcVar psir_x = dataFile.addVar ("psi_r",   ncDouble, psi_d);
      psir_x.putVar (Psi_z);
      NcVar pspr_x = dataFile.addVar ("psi_pr",  ncDouble, psi_d);
      pspr_x.putVar (Psp_z);
      NcVar pppr_x = dataFile.addVar ("psi_ppr", ncDouble, psi_d);
      pppr_x.putVar (Ppp_z);
      NcVar curr_x = dataFile.addVar ("j_r",     ncDouble, psi_d);
      curr_x.putVar (Cur_z);
    }
      catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/RegularTear/RegularTear.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] Psi_y; delete[] Cur_y; delete[] Psi_z; delete[] Cur_z; delete[] Psp_y; delete[] Psp_z;
  delete[] Ppp_y; delete[] Ppp_z;
}

// #########################################
// Function to return equilibrium quantities
// #########################################
void RegularTear::GetEquilibrium (double r, double& q, double& s, double& J, double& Jp, double& lambda)
{
  if (r > 1.)
    {
      J      = 0.;
      Jp     = 0.;
      lambda = 0.;
    }
  else
    {
      double r2  = r*r;
      double r4  = r2*r2;
      double or2 = 1. - r2;
      double nu1 = nu - 1.;
      double nu2 = nu - 2.;
      
      double f  = (1. - pow (or2, nu)) /q0/nu;
      double fp = 2. * r * pow (or2, nu1) /q0;
      
      if (r < 0.01)
	{
	  q = q0 * (1. + 0.5*nu1*r2 + (0.25*nu1*nu1 - nu1*nu2/6.)*r4);
	  s = q0 * (nu1*r2 + 4.*(0.25*nu1*nu1 - nu1*nu2/6.)*r4) /q;
	}
      else
	{
	  q = r2 /f;
	  s = 2. - r*fp /f;
	}
      
      J      = (2./q0) * pow (or2, nu1);
      Jp     = - (4./q0) * nu1 * r * pow (or2, nu2);
      lambda = - q * Jp /s;
    }
}

// #########################################
// Function to return value of safety-factor
// #########################################
double RegularTear::Getq (double r)
{
  double r2  = r*r;
  double r4  = r2*r2;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;
  double q;

  double f  = (1. - pow (or2, nu)) /q0/nu;
  double fp = 2. * r * pow (or2, nu1) /q0;
  
  if (r < 0.01)
    {
      q = q0 * (1. + 0.5*nu1*r2 + (0.25*nu1*nu1 - nu1*nu2/6.)*r4);
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
double RegularTear::Getqp (double r)
{
  double r2  = r*r;
  double r3  = r*r2;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;
  double qp;

  double f  = (1. - pow (or2, nu)) /q0/nu;
  double fp = 2. * r * pow (or2, nu1) /q0;
  
  if (r < 0.01)
    {
      qp = q0 * (nu1*r + 4.*(0.25*nu1*nu1 - nu1*nu2/6.)*r3);
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
double RegularTear::FindRationalSurface ()
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
	      printf ("RegularTear:FindRationalSurface - Error: rn is NaN\n");
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
double RegularTear::GetDelta (int isurf)
{
  // --------------------------------------------------------------------
  // Launch solution from magnetic axis and integrate to rational surface
  // --------------------------------------------------------------------
  double r, h, t_err;
  int    rept;
  double y[2], err[2], dydr[2];
  
  r     = eps;
  y[0]  = pow (r, mpol);
  y[1]  = mpol * pow (r, mpol-1.);
  h     = h0;
  count = 0;

  printf ("\nLaunching solution from magnetic axis:   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  double q, s, J, Jp, lambda;
  GetEquilibrium (r, q, s, J, Jp, lambda);
  CashKarp45Rhs  (r, y, dydr);
  Psi  (0, isurf) = y[0];
  Psip (0, isurf) = y[1];
  Psipp(0, isurf) = dydr[1];
  Cur  (0, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

  int isave;
  for (int i = 1; i < Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);

      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      Psi  (i, isurf) = y[0];
      Psip (i, isurf) = y[1];
      Psipp(i, isurf) = dydr[1];
      Cur  (i, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

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
  GetEquilibrium (rs, q, s, J, Jp, lambda);
  
  double a11 = 1. - lambda * del * (log (del) - 1.);
  double a12 =  - del;
  double a21 = lambda * log (del);
  double a22 = 1.;
  double jac = a11 * a22 - a12 * a21;

  double Clm = (  a22 * y[0] - a12 * y[1]) /jac;
  double Csm = (- a21 * y[0] + a11 * y[1]) /jac;

  Norm (isurf) = Clm;
  for (int i = 0; i <= isave; i++)
    {
      Psi  (i, isurf) /= Clm;
      Psip (i, isurf) /= Clm;
      Psipp(i, isurf) /= Clm;
      Cur  (i, isurf) /= Clm;
    }
  
  printf ("Stopping solution at rational surface:   r = %11.4e y[0] = %11.4e y[1] = %11.4e Cl = %11.4e Cs = %11.4e\n",
	  r, y[0], y[1], Clm, Csm);

  // ----------------------------------------------------------------------
  // Launch solution from plasma boundary and integrate to rational surface
  // ----------------------------------------------------------------------
  r     = 1.;
  y[0]  = 1.;
  y[1]  = - mpol;
  h     = - h0;
  count = 0;

  GetEquilibrium (r, q, s, J, Jp, lambda);
  CashKarp45Rhs  (r, y, dydr);
  Psi  (Nr-1, isurf) = y[0];
  Psip (Nr-1, isurf) = y[1];
  Psipp(Nr-1, isurf) = dydr[1];
  Cur  (Nr-1, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

  if (Fixed)
    {
      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      y[0]               = 0.;
      y[1]               = - 1.;
      Psi  (Nr-1, isurf) = y[0];
      Psip (Nr-1, isurf) = y[1];
      Psipp(Nr-1, isurf) = dydr[1];
      Cur  (Nr-1, isurf) = Jp * y[0] /(1./q - 1./qs) /r;
    }

  printf ("Launching solution from plasma boundary: r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  for (int i = Nr-2; i > 0; i--)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);

      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      Psi  (i, isurf) = y[0];
      Psip( i, isurf) = y[1];
      Psipp(i, isurf) = dydr[1];
      Cur  (i, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

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

  Norp(isurf) = Clp;
  for (int i = isave; i < Nr; i++)
    {
      Psi  (i, isurf) /= Clp;
      Psip (i, isurf) /= Clp;
      Psipp(i, isurf) /= Clp;
      Cur  (i, isurf) /= Clp;
    }

  printf ("Stopping solution at rational surface:   r = %11.4e y[0] = %11.4e y[1] = %11.4e Cl = %11.4e Cs = %11.4e\n",
	  r, y[0], y[1], Clp, Csp);

  // ---------------------------------
  // Calculate tearing stability index
  // ---------------------------------
  double Delta_ = rs * (Csp /Clp - Csm /Clm);

  return Delta_;
}

// #########################################################
// Function to calculate regularized tearing stability index
// #########################################################
double RegularTear::GetRegularDelta (int isurf)
{
  // --------------------------------------------------------------------
  // Launch solution from magnetic axis and integrate to rational surface
  // --------------------------------------------------------------------
  double r, h, t_err;
  int    rept;
  double y[2], err[2], dydr[2];
  
  r     = eps;
  y[0]  = pow (r, mpol)           /Norm(isurf);
  y[1]  = mpol * pow (r, mpol-1.) /Norm(isurf);
  h     = h0;
  count = 0;

  printf ("\nLaunching solution from magnetic axis:   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  double q, s, J, Jp, lambda;
  GetEquilibrium (r, q, s, J, Jp, lambda);
  CashKarp45Rhs  (r, y, dydr);
  Psir  (0, isurf) = y[0];
  Psipr (0, isurf) = y[1];
  Psippr(0, isurf) = dydr[1];
  Curr  (0, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

  int isave;
  for (int i = 1; i < Nr; i++)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r < rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);

      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      Psir  (i, isurf) = y[0];
      Psipr (i, isurf) = y[1];
      Psippr(i, isurf) = dydr[1];
      Curr  (i, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

      if (rr[i+1] > rs - 3.*sigma)
	{
	  isave = i;
	  break;
	}
    }

  do
    {
      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r < rs - 3.*sigma);
  CashKarp45Fixed (2, r, y, err, rs - 3.*sigma - r);

  double y0m = y[0];
  double y1m = y[1];

  CashKarp45Rhs (r, y, dydr);
  double y2m = dydr[1];
  
  printf ("Stopping solution at rational surface:   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  // ----------------------------------------------------------------------
  // Launch solution from plasma boundary and integrate to rational surface
  // ----------------------------------------------------------------------
  r     = 1.;
  y[0]  = 1.     /Norp(isurf);
  y[1]  = - mpol /Norp(isurf);
  h     = - h0;
  count = 0;

  GetEquilibrium (r, q, s, J, Jp, lambda);
  CashKarp45Rhs  (r, y, dydr);
  Psir  (Nr-1, isurf) = y[0];
  Psipr (Nr-1, isurf) = y[1];
  Psippr(Nr-1, isurf) = dydr[1];
  Curr  (Nr-1, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

  if (Fixed)
    {
      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      y[0] = 0.;
      y[1] = - 1. /Norp(isurf);

      Psir  (Nr-1, isurf) = y[0];
      Psipr (Nr-1, isurf) = y[1];
      Psippr(Nr-1, isurf) = dydr[1];
      Curr  (Nr-1, isurf) = Jp * y[0] /(1./q - 1./qs) /r;
    }

  printf ("Launching solution from plasma boundary: r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  for (int i = Nr-2; i > 0; i--)
    {
      do
	{
	  CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (r > rr[i]);
      CashKarp45Fixed (2, r, y, err, rr[i] - r);

      GetEquilibrium (r, q, s, J, Jp, lambda);
      CashKarp45Rhs  (r, y, dydr);
      Psir  (i, isurf) = y[0];
      Psipr (i, isurf) = y[1];
      Psippr(i, isurf) = dydr[1];
      Curr  (i, isurf) = Jp * y[0] /(1./q - 1./qs) /r;

      if (rr[i-1] < rs + 3.*sigma)
	{
	  isave = i;
	  break;
	}
    }
  
  do
    {
      CashKarp45Adaptive (2, r, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (r > rs + 3.*sigma);
  CashKarp45Fixed (2, r, y, err, rs + 3.*sigma - r);

  double y0p = y[0];
  double y1p = y[1];

  CashKarp45Rhs (r, y, dydr);
  double y2p = dydr[1];

  printf ("Stopping solution at rational surface:   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

  // ---------------------------------
  // Calculate tearing stability index
  // ---------------------------------
  double C0 = (5./3.) * sigma*sigma * (y2p - y2m);
  double C3 = 0.5 * sigma*sigma * (y2p + y2m);
  double C4 = 0.5 * sigma * (y1p + y1m) - 0.5 * log (10.) * C0;
  double C2 = sigma * (y1p - y1m) - 6. * C3;
  double C5 = 0.5 * (y0p + y0m) - 1.5 * C2 - 4.5 * C3;
  double C1 = (y0p - y0m) - 6. * C4 - 2. * C0 * (1.5 * (log (10.) - 2.) + atan (3.));

  printf ("Fit                                  :   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  rs - 3.*sigma, GetPsiFit (-3., C0, C1, C2, C3, C4, C5),  GetPsipFit (-3., C0, C1, C2, C3, C4) /sigma);
  printf ("Fit                                  :   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  rs + 3.*sigma, GetPsiFit (+3., C0, C1, C2, C3, C4, C5),  GetPsipFit (+3., C0, C1, C2, C3, C4) /sigma);
  
  printf ("\nC0 = %11.4e C1 = %11.4e C2 = %11.4e C3 = %11.4e C4 = %11.4e C5 = %11.4e\n", C0, C1, C2, C3, C4, C5);

  for (int i = 0; i < Nr; i++)
    {
      if (rr[i] > rs - 3.*sigma && rr[i] < rs + 3.*sigma)
	{
	  double X = (rr[i] - rs) /sigma;

	  Psir  (i, isurf) = GetPsiFit  (X, C0, C1, C2, C3, C4, C5);  
	  Psipr (i, isurf) = GetPsipFit (X, C0, C1, C2, C3, C4)     /sigma;
	  Psippr(i, isurf) = GetPsippFit(X, C0, C1, C2, C3)         /sigma/sigma;  
	  Curr  (i, isurf) = GetCurFit  (X, C0, C1, C2, C3, C4, C5);  
	}
    }
  
  double Delta_ = rs * C2 /C5 /sigma;

  return Delta_;
}

// ##############################
// Function to return fit for Psi
// ##############################
double RegularTear::GetPsiFit (double X, double C0, double C1, double C2, double C3, double C4, double C5)
{
  return  C0 * (0.5 * X * (log (1. + X*X) - 2.) + atan (X))
        + 0.5 * C1 * gsl_sf_erf (X /sqrt(2.))
        + C2 * (0.5 * X * gsl_sf_erf (X /sqrt(2.)) + exp (-X*X/2.) /sqrt (2.*M_PI))
        + 0.5 * C3 * X*X + C4 * X  + C5;
}

// ###############################
// Function to return fit for Psip
// ###############################
double RegularTear::GetPsipFit (double X, double C0, double C1, double C2, double C3, double C4)
{
  return 0.5 * C0 * log (1. + X*X) + C1 * exp (-X*X/2.) /sqrt (2.*M_PI) + 0.5 * C2 * gsl_sf_erf (X /sqrt(2.)) + C3 * X + C4;
}

// ################################
// Function to return fit for Psipp
// ################################
double RegularTear::GetPsippFit (double X, double C0, double C1, double C2, double C3)
{
  return C0 * X /(1. + X*X) - C1 * X * exp(-X*X/2.) /sqrt (2.*M_PI) + C2 * exp (-X*X/2.) /sqrt (2.*M_PI) + C3;
}

// ##################################
// Function to return fit for Current
// ##################################
double RegularTear::GetCurFit (double X, double C0, double C1, double C2, double C3, double C4, double C5)
{
  double psi   = GetPsiFit   (X, C0, C1, C2, C3, C4, C5);
  double psip  = GetPsipFit  (X, C0, C1, C2, C3, C4);
  double psipp = GetPsippFit (X, C0, C1, C2, C3);

  double r = rs + sigma * X;
  
  return psipp /sigma/sigma + psip /r /sigma - mpol*mpol * psi /r/r;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void RegularTear::CashKarp45Rhs (double r, double* y, double* dydr)
{
  // y[0] = psi
  // y[1] = psi'
  
  double q, s, J, Jp, lambda;
  
  GetEquilibrium (r, q, s, J, Jp, lambda);
  
  dydr[0] = y[1];
  dydr[1] = - y[1] /r + mpol*mpol * y[0] /r/r + Jp * y[0] /r /(1./q - 1./qs);
}

