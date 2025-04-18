// Tear.cpp

#include "Tear.h"

// ###########
// Constructor
// ###########
Tear::Tear ()
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
  string JSONFilename = "../Inputs/Tear.json";
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

  // ------------
  // Sanity check
  // ------------
  if (q0 < 0.)
    {
      printf ("Tear:: Error - q0 cannot be negative\n");
      exit (1);
    }
  if (qa < 2.*q0)
    {
      printf ("Tear:: Error - qa cannot be less than 2*q0\n");
      exit (1);
    }
 if (NTOR < 1)
    {
      printf ("Tear:: Error - NTOR must be positive\n");
      exit (1);
    }
  if (eps <= 0.)
    {
      printf ("Tear:: Error - eps must be positive\n");
      exit (1);
    }
  if (del <= 0.)
    {
      printf ("Tear:: Error - del must be positive\n");
      exit (1);
    }
  if (Nr < 2)
    {
      printf ("Tear:: Error - Nr cannot be less than two\n");
      exit (1);
    }
  if (acc <= 0.)
    {
      printf ("Tear:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("Tear:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("Tear:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("Tear:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("Tear:: Error - hmax must exceed hmin\n");
      exit (1);
    }

  ntor = double (NTOR);
  nu   = qa /q0;

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  
  printf ("\nClass TEAR::\n\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ntor = %-2d    q0  = %10.3e qa  = %10.3e Fixed = %1d\n",
	  NTOR, q0, qa, Fixed);
  printf ("Nr   = %4d eps  = %10.3e del = %10.3e acc   = %10.3e h0 = %10.3e hmax = %10.3e\n",
	  Nr, eps, del, acc, h0, hmax);
}

// ##########
// Destructor
// ##########
Tear::~Tear ()
{
}

// #########################
// Function to solve problem
// #########################
void Tear::Solve ()
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
  mres  = new int[nres];
  rres  = new double[nres];
  Dres  = new double[nres];
  Psi.resize (Nr, nres);
    
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
 
      if (rs > 0. && rs < 1.)
	{
	  printf ("\nFound rational surface: mpol = %-2d ntor = %-2d rs = %11.4e residual = %11.4e\n",
		  mmin + isurf, NTOR, rs, fabs (Getq(rs) - qs));
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

      // ----------
      // Store data
      // ----------
      mres[isurf] = mmin + isurf;
      rres[isurf] = rs;
      Dres[isurf] = Delta;
    }
 
  // -----------------
  // Write netcdf file
  // -----------------
  WriteNetcdf ();
  
  // -------
  // Cleanup
  // -------
  delete[] rr;   delete[] ss;   delete[] JJ;   delete[] JJp; delete[] lvals;
  delete[] mres; delete[] rres; delete[] Dres;
}

// #####################################
// Function to write data to netcdf file
// #####################################
void Tear::WriteNetcdf ()
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
      NcFile dataFile ("../Outputs/Tear/Tear.nc", NcFile::replace);
       
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
      NcVar D_x   = dataFile.addVar ("Dres",   ncDouble, s_d);
      D_x.putVar (Dres);
      NcVar psi_x = dataFile.addVar ("psi",    ncDouble, psi_d);
      psi_x.putVar (Psi_y);
    }
      catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Tear/Tear.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] Psi_y;
}

// #########################################
// Function to return equilibrium quantities
// #########################################
void Tear::GetEquilibrium (double r, double& q, double& s, double& J, double& Jp, double& lambda)
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

// #########################################
// Function to return value of safety-factor
// #########################################
double Tear::Getq (double r)
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
double Tear::Getqp (double r)
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
double Tear::FindRationalSurface ()
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
	      printf ("Tear:FindRationalSurface - Error: rn is NaN\n");
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
double Tear::GetDelta (int isurf)
{
  // --------------------------------------------------------------------
  // Launch solution from magnetic axis and integrate to rational surface
  // --------------------------------------------------------------------
  double r, h, t_err;
  int    rept;
  double y[2], err[2];
  
  r     = eps;
  y[0]  = pow (r, mpol);
  y[1]  = mpol * pow (r, mpol-1.);
  h     = h0;
  count = 0;

  printf ("\nLaunching solution from magnetic axis:   r = %11.4e y[0] = %11.4e y[1] = %11.4e\n",
	  r, y[0], y[1]);

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

  Psi(Nr-1, isurf) = 1.;

  if (Fixed)
    {
      y[0]             = 0.;
      y[1]             = - 1.;
      Psi(Nr-1, isurf) = 0.;
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

  printf ("Stopping solution at rational surface:   r = %11.4e y[0] = %11.4e y[1] = %11.4e Cl = %11.4e Cs = %11.4e\n",
	  r, y[0], y[1], Clp, Csp);

  // ---------------------------------
  // Calculate tearing stability index
  // ---------------------------------
  double Delta_ = rs * (Csp /Clp - Csm /Clm);

  return Delta_;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Tear::CashKarp45Rhs (double r, double* y, double* dydr)
{
  // y[0] = psi
  // y[1] = dpsi/dr
  
  double q, s, J, Jp, lambda;

  GetEquilibrium (r, q, s, J, Jp, lambda);

  dydr[0] = y[1];
  dydr[1] = - y[1] /r + mpol*mpol * y[0] /r/r + Jp * y[0] /r /(1./q - 1./qs); 
 }

