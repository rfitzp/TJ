// TearX.cpp

#include "TearX.h"

// ###########
// Constructor
// ###########
TearX::TearX ()
{
  // --------------------------------------------
  // Ensure that directory ../Outputs/TearX exits
  // --------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/TearX"))
    {
      exit (1);
    }

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/TearX.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  NTOR    = JSONData["NTOR"]   .get<int>    ();
  q0      = JSONData["q0"]     .get<double> ();
  nu      = JSONData["nu"]     .get<double> ();
  alpha   = JSONData["alpha"]  .get<double> ();
  Delta_q = JSONData["Delta_q"].get<double> ();
  r_q     = JSONData["r_q"]    .get<double> ();
  delta_q = JSONData["delta_q"].get<double> ();

  B0   = JSONData["B0"]  .get<double> ();
  R0   = JSONData["R0"]  .get<double> ();
  a    = JSONData["a"]   .get<double> ();
  M    = JSONData["M"]   .get<double> ();
  chiE = JSONData["chiE"].get<double> ();
  chip = JSONData["chip"].get<double> ();
  Z    = JSONData["Z"]   .get<double> ();
 
  ne_0     = JSONData["ne_0"]    .get<double> ();
  ne_1     = JSONData["ne_1"]    .get<double> ();
  Delta_ne = JSONData["Delta_ne"].get<double> ();
  nu_ne    = JSONData["nu_ne"]   .get<double> ();
  r_ne     = JSONData["r_ne"]    .get<double> ();
  delta_ne = JSONData["delta_ne"].get<double> ();

  Te_0     = JSONData["Te_0"]    .get<double> ();
  Te_1     = JSONData["Te_1"]    .get<double> ();
  Delta_Te = JSONData["Delta_Te"].get<double> ();
  nu_Te    = JSONData["nu_Te"]   .get<double> ();
  r_Te     = JSONData["r_Te"]    .get<double> ();
  delta_Te = JSONData["delta_Te"].get<double> ();

  Ti_0     = JSONData["Ti_0"]    .get<double> ();
  Ti_1     = JSONData["Ti_1"]    .get<double> ();
  Delta_Ti = JSONData["Delta_Ti"].get<double> ();
  nu_Ti    = JSONData["nu_Ti"]   .get<double> ();
  r_Ti     = JSONData["r_Ti"]    .get<double> ();
  delta_Ti = JSONData["delta_Ti"].get<double> ();

  alphaE   = JSONData["alphaE"]  .get<double> ();
  omegaE_0 = JSONData["omegaE_0"].get<double> ();
  nu_E     = JSONData["nu_E"]    .get<double> ();
  
  eps    = JSONData["eps"]   .get<double> ();
  del    = JSONData["del"]   .get<double> ();
  EPS    = JSONData["EPS"]   .get<double> ();
  Psimax = JSONData["Psimax"].get<double> ();
  Nr     = JSONData["Nr"]    .get<int>    ();

  nu_sta = JSONData["nu_sta"].get<double> ();
  nu_end = JSONData["nu_end"].get<double> ();
  nu_num = JSONData["nu_num"].get<int>    ();
  m_min  = JSONData["m_min"] .get<int>    ();
  r_min  = JSONData["r_min"] .get<double> ();
  c_min  = JSONData["c_min"] .get<double> ();
  
  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

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
  if (nu_num < 2)
    {
      printf ("TearX:: Error - nu_num must be greater than 1\n");
      exit (1);
    }
  if (nu_ne < 1.)
    {
      printf ("TearX:: Error - nu_ne must be greater than 1\n");
      exit (1);
    }
  if (nu_Te < 1.)
    {
      printf ("TearX:: Error - nu_Te must be greater than 1\n");
      exit (1);
    }
  if (nu_Ti < 1.)
    {
      printf ("TearX:: Error - nu_Ti must be greater than 1\n");
      exit (1);
    }
 if (nu_E < 1.)
    {
      printf ("TearX:: Error - nu_E must be greater than 1\n");
      exit (1);
    }

  ntor = double (NTOR);
  qc   = nu * q0;

  e         = GSL_CONST_MKSA_ELECTRON_CHARGE;
  m_e       = GSL_CONST_MKSA_MASS_ELECTRON;
  m_p       = GSL_CONST_MKSA_MASS_PROTON;
  epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
  mu_0      = GSL_CONST_MKSA_VACUUM_PERMEABILITY;

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\nClass TEARX::\n\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("ntor   = %-2d         q0       = %-10.3e nu       = %-10.3e qc    = %-10.3e alpha  = %-10.3e Delta_q  = %-10.3e r_q = %-10.3e delta_q = %-10.3e\n",
	  NTOR, q0, nu, qc, alpha, Delta_q, r_q, delta_q);
  printf ("B0     = %-10.3e R0       = %-10.3e a        = %-10.3e M     = %-10.3e chiE   = %-10.3e chip     = %-10.3e Z   = %-10.3e\n",
	  B0, R0, a, M, chiE, chip, Z);
  printf ("ne_0   = %-10.3e ne_1     = %-10.3e Delta_ne = %-10.3e nu_ne = %-10.3e r_ne   = %-10.3e delta_ne = %-10.3e\n",
	  ne_0, ne_1, Delta_ne, nu_ne, r_ne, delta_ne);
  printf ("Te_0   = %-10.3e Te_1     = %-10.3e Delta_Te = %-10.3e nu_Te = %-10.3e r_Te   = %-10.3e delta_Te = %-10.3e\n",
	  Te_0, Te_1, Delta_Te, nu_Te, r_Te, delta_Te);
  printf ("Ti_0   = %-10.3e Ti_1     = %-10.3e Delta_Ti = %-10.3e nu_Ti = %-10.3e r_Ti   = %-10.3e delta_Ti = %-10.3e\n",
	  Ti_0, Ti_1, Delta_Ti, nu_Ti, r_Ti, delta_Ti);
  printf ("alphaE = %-10.3e omegaE_0 = %-10.3e nu_E     = %-10.3e\n",
	  alphaE, omegaE_0, nu_E);
  printf ("\nNr     = %-5d      eps      = %-10.3e del      = %-10.3e EPS   = %-10.3e Psimax = %-10.3e\n",
	  Nr, eps, del, EPS, Psimax);
  printf ("nu_sta = %-10.3e nu_end   = %-10.3e nu_num   = %-5d      m_min = %-3d        r_min  = %-10.3e c_min    = %-10.3e\n",
	  nu_sta, nu_end, nu_num, m_min, r_min, c_min);
  printf ("acc    = %-10.3e h0       = %-10.3e hmax     = %-10.3e\n",
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
// Function to scan nu
// ###################
void TearX::Scannu ()
{
  FILE* file = OpenFilew ("../Outputs/TearX/Scannu.txt");
  fclose (file);

  for (int i = 0; i < nu_num; i++)
    {
      double _nu = nu_sta + (nu_end - nu_sta) * double (i) /double (nu_num - 1);

      Setnu (_nu);

      Solve (1);
    }
}

// #########################
// Function to solve problem
// #########################
void TearX::Solve (int flg)
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
  Q_spline = gsl_spline_alloc (gsl_interp_cspline, Nr);

  q_acc = gsl_interp_accel_alloc ();
  P_acc = gsl_interp_accel_alloc ();
  R_acc = gsl_interp_accel_alloc ();
  Q_acc = gsl_interp_accel_alloc ();

  gsl_spline_init (q_spline, rr,  qq,  Nr);
  gsl_spline_init (P_spline, PSI, rr,  Nr);
  gsl_spline_init (R_spline, rr,  rho, Nr);
  gsl_spline_init (Q_spline, rr,  PSI, Nr);

  r95 = gsl_spline_eval (P_spline, 0.95, P_acc);
  R95 = gsl_spline_eval (R_spline, r95,  R_acc);
  q95 = Getq (r95);
  s95 = GetShear (r95);
  S95 = R95 * gsl_spline_eval_deriv (q_spline, r95, q_acc) /q95 /gsl_spline_eval_deriv (R_spline, r95, R_acc);

  rmax = gsl_spline_eval (P_spline, Psimax, P_acc);
  qmax = gsl_spline_eval (q_spline, rmax,   P_acc);

  if (flg == 0)
    printf ("\nr95  = %-10.3e R95  = %-10.3e q95  = %-10.3e s95 = %-10.3e S95 = %10.3e\n",
	    r95, R95, q95, s95, S95);
  if (flg == 0)
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
  Drres = new double[nres];
  Dires = new double[nres];
  Idres = new double[nres];
  Psi.resize (Nr, nres);
    
  ss    = new double[Nr];
  JJ    = new double[Nr];
  JJp   = new double[Nr];
  JJbs  = new double[Nr];
  lvals = new double[Nr];

  ne    = new double[Nr];
  dnedr = new double[Nr];
  Te    = new double[Nr];
  dTedr = new double[Nr];
  Ti    = new double[Nr];
  dTidr = new double[Nr];

  waste = new double[Nr];
  wasti = new double[Nr];
  wE    = new double[Nr];

  tauR  = new double[Nr];
  tauA  = new double[Nr];
  tauE  = new double[Nr];
  taup  = new double[Nr];
  S     = new double[Nr];
  cbeta = new double[Nr];
  dbeta = new double[Nr];

  QE    = new double[Nr];
  Qe    = new double[Nr];
  Qi    = new double[Nr];
  D     = new double[Nr];
  P_E   = new double[Nr];
  P_phi = new double[Nr];
  Scale = new double[Nr];
  
  // --------------------------------
  // Calculate equilibrium quantities
  // --------------------------------
  for (int i = 0; i < Nr; i++)
    {
      double q, s, J, Jp, lambda;
      
      GetEquilibrium (rr[i], q, s, J, Jp, lambda);

      ss [i] = s;
      JJ [i] = J;
      JJp[i] = Jp;

      ne   [i] = Getne    (rr[i]);
      dnedr[i] = Getdnedr (rr[i]);
      Te   [i] = GetTe    (rr[i]);
      dTedr[i] = GetdTedr (rr[i]);
      Ti   [i] = GetTi    (rr[i]);
      dTidr[i] = GetdTidr (rr[i]);
      JJbs [i] = GetJbs   (rr[i]);

      waste[i] = Getwaste (rr[i]);
      wasti[i] = Getwasti (rr[i]);
      wE   [i] = GetwE    (rr[i]);

      tauR [i] = GettauR (rr[i]);
      tauA [i] = GettauA (rr[i]);
      tauE [i] = GettauE (rr[i]);
      taup [i] = Gettaup (rr[i]);
      S    [i] = tauR[i] /tauA[i];
      cbeta[i] = Getcbeta (rr[i]);
      dbeta[i] = Getdbeta (rr[i]);

      QE   [i] = GetQE    (rr[i]);
      Qe   [i] = GetQe    (rr[i]);
      Qi   [i] = GetQi    (rr[i]);
      D    [i] = GetD     (rr[i]);
      P_E  [i] = GetP_E   (rr[i]);
      P_phi[i] = GetP_phi (rr[i]);
      Scale[i] = GetScale (rr[i]);

      if (i == 0)
	lvals[i] = 4. * q0 * (nu - 1.) /(q0 * (nu - 1.) + 2. * alpha);
      else
	lvals[i] = rr[i]*lambda;
    }
  waste[0] = waste[1];
  wasti[0] = wasti[1];
  wE   [0] = wE   [1];
  JJbs [0] = JJbs [1];

  double smin = 1.e6, rmin;
  for (int i = 0; i < Nr; i++)
    if (ss[i] < smin && rr[i] > r_min)
      {
	smin = ss[i];
	rmin = rr[i];
      }
	
  if (flg == 0)
    printf ("s_min = %10.3e r_min = %-10.3e ne(1) = %-10.3e Te(1) = %-10.3e Ti(1) = %-10.3e\n\n",
	    smin, rmin, ne[Nr-1], Te[Nr-1], Ti[Nr-1]);

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
      
      // ------------------------------
      // Calculate layer response index
      // ------------------------------
      double QEs = GetQE    (rs);
      double Qes = GetQe    (rs);
      double Qis = GetQi    (rs);
      double cbs = Getcbeta (rs);
      double Ds  = GetD     (rs);
      double PEs = GetP_E   (rs);
      double Pps = GetP_phi (rs);
      double Sc  = GetScale (rs);
      double Delta_r, Delta_i, pmax, Ideal;

      if (rs > r_min && mpol >= m_min)
	{
	  FourField fourfield (1);

	  if (cbs < c_min)
	    fourfield.GetDelta3 (QEs, Qes, Qis, cbs, Ds, PEs, Pps, Delta_r, Delta_i, pmax);
	  else
	    fourfield.GetDelta4 (QEs, Qes, Qis, cbs, Ds, PEs, Pps, Delta_r, Delta_i, pmax);
	  if (flg != 0)
	    printf ("QE = %9.2e Qe = %9.2e Qi = %9.2e cbeta = %9.2e D = %9.2e PE = %9.2e Pp = %9.2e Delta = (%9.2e, %9.2e) pmax = %9.2e\n",
		    QEs, Qes, Qis, cbs, Ds, PEs, Pps, Delta_r, Delta_i, pmax);
	  
	  Delta_r *= Sc;
	  Delta_i *= Sc;
	  
	  Ideal = sqrt (Delta_r*Delta_r + Delta_i*Delta_i) /fabs (Delta);
	  
	  if (flg == 0)
	    printf ("mpol = %-2d ntor = %-2d rs = %-10.3e s = %-10.3e E_ss = %10.3e Delta = (%10.3e, %10.3e) Ideal = %10.3e\n",
		    mmin + isurf, NTOR, rs, sr, Delta, Delta_r, Delta_i, Ideal);
	}
      else
	{
	  Delta_r = 0.;
	  Delta_i = 0.;
	  Ideal   = -1.;
	}

      // ----------
      // Store data
      // ----------
      mres [isurf] = mmin + isurf;
      rres [isurf] = rs;
      Pres [isurf] = gsl_spline_eval (Q_spline, rres[isurf], Q_acc);
      sres [isurf] = sr;
      Dres [isurf] = Delta;
      Drres[isurf] = Delta_r;
      Dires[isurf] = Delta_i;
      Idres[isurf] = Ideal;
    }

  // -----------------
  // Write netcdf file
  // -----------------
  if (flg == 0)
    WriteNetcdf ();

  // ---------------
  // Write scan data
  // ---------------
  if (flg != 0)
    {
      FILE* file = OpenFilea ("../Outputs/TearX/Scannu.txt");
      
      printf ("nu = %-10.3e qa/qc = %-10.3e q95 = %-10.3e nres = %-3d m = (%d, %d)\n",
	      nu, qq[Nr-1]/qc, q95, nres, mres[0], mres[nres-1]);
	      
      for (int i = 0; i < nres; i++)
	if (rres[i] > r_min && mres[i] >= m_min)
	  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e %3d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
		   q0, nu, qc, r95, q95, mres[i], rres[i], sres[i], Dres[i], Pres[i], Drres[i], Dires[i], Idres[i]);
      fclose (file);
    }
  
  // -------
  // Cleanup
  // -------
  delete[] rr;   delete[] ss;   delete[] JJ;    delete[] JJp;   delete[] lvals;
  delete[] mres; delete[] rres; delete[] Dres;  delete[] PSI;   delete[] sres;
  delete[] rho;  delete[] Pres; delete[] Drres; delete[] Dires; delete[] Idres;

  delete[] ne;    delete[] dnedr; delete[] Te;    delete[] dTedr; delete[] Ti;   delete[] dTidr;
  delete[] waste; delete[] wasti; delete[] wE;    delete[] tauR;  delete[] tauA; delete[] tauE;
  delete[] taup;  delete[] S;     delete[] cbeta; delete[] dbeta; delete[] JJbs; delete[] Scale;
  delete[] QE;    delete[] Qe;    delete[] Qi;    delete[] D;     delete[] P_E;  delete[] P_phi; 

  gsl_spline_free (q_spline);    gsl_spline_free (P_spline);    gsl_spline_free (R_spline);    gsl_spline_free (Q_spline);
  gsl_interp_accel_free (q_acc); gsl_interp_accel_free (P_acc); gsl_interp_accel_free (R_acc); gsl_interp_accel_free (Q_acc);
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
      NcVar ne_x  = dataFile.addVar ("ne",     ncDouble, r_d);
      ne_x.putVar (ne);
      NcVar nep_x = dataFile.addVar ("dnedr",  ncDouble, r_d);
      nep_x.putVar (dnedr);
      NcVar Te_x  = dataFile.addVar ("Te",     ncDouble, r_d);
      Te_x.putVar (Te);
      NcVar Tep_x = dataFile.addVar ("dTedr",  ncDouble, r_d);
      Tep_x.putVar (dTedr);
      NcVar Ti_x  = dataFile.addVar ("Ti",     ncDouble, r_d);
      Ti_x.putVar (Ti);
      NcVar Tip_x = dataFile.addVar ("dTidr",  ncDouble, r_d);
      Tip_x.putVar (dTidr);
      NcVar jbs_x = dataFile.addVar ("Jbs",   ncDouble, r_d);
      jbs_x.putVar (JJbs);

      NcVar we_x = dataFile.addVar  ("omegae", ncDouble, r_d);
      we_x.putVar (waste);
      NcVar wi_x = dataFile.addVar  ("omegai", ncDouble, r_d);
      wi_x.putVar (wasti);
      NcVar wE_x = dataFile.addVar  ("omegaE", ncDouble, r_d);
      wE_x.putVar (wE);

      NcVar tr_x = dataFile.addVar  ("tauR",   ncDouble, r_d);
      tr_x.putVar (tauR);
      NcVar ta_x = dataFile.addVar  ("tauA",   ncDouble, r_d);
      ta_x.putVar (tauA);
      NcVar te_x = dataFile.addVar  ("tauE",   ncDouble, r_d);
      te_x.putVar (tauE);
      NcVar tp_x = dataFile.addVar  ("taup",   ncDouble, r_d);
      tp_x.putVar (taup);
      NcVar S_x = dataFile.addVar   ("S",      ncDouble, r_d);
      S_x.putVar (S);
      NcVar cb_x = dataFile.addVar  ("c_beta", ncDouble, r_d);
      cb_x.putVar (cbeta);
      NcVar db_x = dataFile.addVar  ("d_beta", ncDouble, r_d);
      db_x.putVar (dbeta);

      NcVar qE_x = dataFile.addVar  ("Q_E",    ncDouble, r_d);
      qE_x.putVar (QE);
      NcVar qe_x = dataFile.addVar  ("Q_e",    ncDouble, r_d);
      qe_x.putVar (Qe);
      NcVar qi_x = dataFile.addVar  ("Q_i",    ncDouble, r_d);
      qi_x.putVar (Qi);
      NcVar d_x = dataFile.addVar   ("D",      ncDouble, r_d);
      d_x.putVar (D);
      NcVar pe_x = dataFile.addVar  ("P_E",    ncDouble, r_d);
      pe_x.putVar (P_E);
      NcVar pp_x = dataFile.addVar  ("P_phi",  ncDouble, r_d);
      pp_x.putVar (P_phi);
      NcVar sc_x = dataFile.addVar  ("Scale",  ncDouble, r_d);
      sc_x.putVar (Scale);
      
      NcVar m_x   = dataFile.addVar ("mres",        ncInt,    s_d);
      m_x.putVar (mres);
      NcVar rs_x  = dataFile.addVar ("rres",        ncDouble, s_d);
      rs_x.putVar (rres);
      NcVar Ps_x  = dataFile.addVar ("Pres",        ncDouble, s_d);
      Ps_x.putVar (Pres);
      NcVar ss_x  = dataFile.addVar ("sres",        ncDouble, s_d);
      ss_x.putVar (sres);
      NcVar D_x   = dataFile.addVar ("Dres",        ncDouble, s_d);
      D_x.putVar (Dres);
      NcVar Dr_x  = dataFile.addVar ("Delta_r_res", ncDouble, s_d);
      Dr_x.putVar (Drres);
      NcVar Di_x  = dataFile.addVar ("Delta_i_res", ncDouble, s_d);
      Di_x.putVar (Dires);
      NcVar id_x  = dataFile.addVar ("Ideal",       ncDouble, s_d);
      id_x.putVar (Idres);

      NcVar psi_x = dataFile.addVar ("psi", ncDouble, psi_d);
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
  q  = Getq (r);
  s  = GetShear (r);

  if (r > 1.)
    {
      J  = 0.;
      Jp = 0.;
    }
  else
    {
      J  = GetJ (r);
      Jp = GetJp (r);
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
  double lg  = log (fabs (1. - r2) + EPS);
 
  double q, f;

  if (r > 1.)
    f = 1. /q0/nu;
  else
    f = (1. - pow (or2, nu)) /q0/nu;
  
  if (r < 0.01)
    {
      q = q0 * (1. + 0.5*nu1*r2 + (0.25*nu1*nu1 - nu1*nu2/6.)*r4) + alpha*r2;
     }
  else
    {
      q = r2 /f - alpha * lg - 0.5 * Delta_q * (1. + tanh ((r - r_q) /delta_q));
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
  double lg  = fabs (1. - r2) + EPS;
  double qp, f, fp;

  if (r > 1.)
    {
      f  = 1. /q0/nu;
      fp = 0.;
    }
  else
    {
      f  = (1. - pow (or2, nu)) /q0/nu;
      fp = 2. * r * pow (or2, nu1) /q0;
    }
  
  if (r < 0.01)
    {
      qp = q0 * (nu1*r + 4.*(0.25*nu1*nu1 - nu1*nu2/6.)*r3) + 2.*alpha*r;
    }
  else
    {
      qp = 2.*r /f - r2 * fp /f/f + 2.*alpha*r /lg - 0.5 * (Delta_q /delta_q) /cosh ((r - r_q) /delta_q) /cosh ((r - r_q) /delta_q);
    }
  
  return qp;
}

// ############################################
// Function to return cylindrical safety-factor
// ############################################
double TearX::Getq_cyl (double r)
{
  double r2  = r*r;
  double r4  = r2*r2;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;
  double q_cyl;

  double f = (1. - pow (or2, nu)) /q0/nu;
  if (r > 1.)
    {
      f = 1. /q0/nu;
    }
  
  if (r < 0.01)
    {
      q_cyl = q0 * (1. + 0.5*nu1*r2 + (0.25*nu1*nu1 - nu1*nu2/6.)*r4);
     }
  else
    {
      q_cyl = r2 /f;
    }

  return q_cyl;
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

// #############################################
// Function to return normalized current density
// #############################################
double TearX::GetJ (double r)
{
  double r2  = r*r;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
 
  return (2./q0) * pow (or2, nu1);
}

// ######################################################
// Function to return normalized current density gradient
// ######################################################
double TearX::GetJp (double r)
{
  double r2  = r*r;
  double or2 = 1. - r2;
  double nu1 = nu - 1.;
  double nu2 = nu - 2.;

  return - (4./q0) * nu1 * r * pow (or2, nu2);
}

// ###################################################
// Function to return value of electron number density
// ###################################################
double TearX::Getne (double r)
{
  return (ne_0 - ne_1) * pow (1 - r*r, nu_ne) + ne_1
    + 0.5 * Delta_ne * (1. - tanh ((r - r_ne) /delta_ne));
}

// ######################################################
// Function to return gradient of electron number density
// ######################################################
double TearX::Getdnedr (double r)
{
  return nu_ne * (ne_0 - ne_1) * (-2.*r) * pow (1 - r*r, nu_ne - 1.)
    - 0.5 * (Delta_ne /delta_ne) /cosh ((r - r_ne) /delta_ne) /cosh ((r - r_ne) /delta_ne);
}

// ################################################
// Function to return value of electron temperature
// ################################################
double TearX::GetTe (double r)
{
  return (Te_0 - Te_1) * pow (1 - r*r, nu_Te) + Te_1
    + 0.5 * Delta_Te * (1. - tanh ((r - r_Te) /delta_Te));
}

// ###################################################
// Function to return gradient of electron temperature
// ###################################################
double TearX::GetdTedr (double r)
{
  return nu_Te * (Te_0 - Te_1) * (-2.*r) * pow (1 - r*r, nu_Te - 1.)
    - 0.5 * (Delta_Te /delta_Te) /cosh ((r - r_Te) /delta_Te) /cosh ((r - r_Te) /delta_Te);
}

// ###########################################
// Function to return value of ion temperature
// ###########################################
double TearX::GetTi (double r)
{
  return (Ti_0 - Ti_1) * pow (1 - r*r, nu_Ti) + Ti_1
    + 0.5 * Delta_Ti * (1. - tanh ((r - r_Ti) /delta_Ti));
}

// ##############################################
// Function to return gradient of ion temperature
// ##############################################
double TearX::GetdTidr (double r)
{
  return nu_Ti * (Ti_0 - Ti_1) * (-2.*r) * pow (1 - r*r, nu_Ti - 1.)
    - 0.5 * (Delta_Ti /delta_Ti) /cosh ((r - r_Ti) /delta_Ti) /cosh ((r - r_Ti) /delta_Ti);
}

// ############################################
// Function to return bootstrap current density
// ############################################
double TearX::GetJbs (double r)
{
  double eps   = (r*a /R0);
  double j0    = mu_0 * R0*R0 * e /a/a /B0/B0;

  double q_cyl = Getq_cyl (r);
  double ne    = Getne (r);
  double nep   = Getdnedr (r);
  double Te    = GetTe (r);
  double Tep   = GetdTedr (r);
  double Ti    = GetTi (r);
  double Tip   = GetdTidr (r);

  return - j0 * (q_cyl * sqrt (eps) /r) * (2.40 * (Te + Ti) * nep + 0.61 * ne * Tep - 0.41 * ne * Tip);
}

// #################################################
// Function to return electron diamagnetic frequency
// #################################################
double TearX::Getwaste (double r)
{
  double q     = Getq    (r);
  double Te    = GetTe   (r);
  double ne    = Getne   (r);
  double dTedr = GetdTedr(r);
  double dnedr = Getdnedr(r);
  
  return (Te * dnedr /ne + dTedr) * (q /B0/a/a) /r;
}

// ############################################
// Function to return ion diamagnetic frequency
// ############################################
double TearX::Getwasti (double r)
{
  double q     = Getq    (r);
  double Ti    = GetTi   (r);
  double ne    = Getne   (r);
  double dTidr = GetdTidr(r);
  double dnedr = Getdnedr(r);
  
  return - (Ti * dnedr /ne + dTidr) * (q /B0/a/a) /r;
}

// ################################
// Function to return ExB frequency
// ################################
double TearX::GetwE (double r)
{
  double waste = Getwaste (r);  

  return omegaE_0 * pow (1. - r*r, nu_E) + alphaE * waste;
}

// #############################################
// Function to return value of Coulomb logarithm
// #############################################
double TearX::GetLambda (double r)
{
  double ne = Getne (r);
  double Te = GetTe (r);

  return 24. - log (sqrt (ne/1.e6) /Te);
}

// ####################################################
// Function to return value of resistive diffusion time
// ####################################################
double TearX::GettauR (double r)
{
  double ne     = Getne (r);
  double Te     = GetTe (r);
  double Lambda = GetLambda (r);

  double tauei = 6. * sqrt (2.) * pow (M_PI, 1.5) * epsilon_0*epsilon_0 * sqrt (m_e) * pow (Te, 1.5)
    /Z /Lambda /pow (e, 2.5) /ne;

  double epsr = r * a /R0;
  double fp   = 1. - 1.46 * sqrt (epsr) + 0.46 * pow (epsr, 1.5);

  double etap = m_e /1.96 /ne /e/e /tauei /fp;

  return mu_0 * r*r * a*a /etap;
}

// #######################################
// Function to return value of Alfven time
// #######################################
double TearX::GettauA (double r)
{
  double ne = Getne (r);
  
  return R0 * sqrt (mu_0 * M * m_p * ne) /B0;
}

// ###################################################
// Function to return value of energy confinement time
// ###################################################
double TearX::GettauE (double r)
{
  return r*r * a*a /chiE;
}

// #####################################################
// Function to return value of momentum confinement time
// #####################################################
double TearX::Gettaup (double r)
{
  return r*r * a*a /chip;
}

// ##################################
// Function to return value of c_beta
// ##################################
double TearX::Getcbeta (double r)
{
  double ne = Getne (r);
  double Te = GetTe (r);
  double Ti = GetTi (r);
  
  double beta = (5./3.) * mu_0 * ne * (Te + Ti) * e /B0/B0;

  return sqrt (beta /(1. + beta));
}

// ##################################
// Function to return value of d_beta
// ##################################
double TearX::Getdbeta (double r)
{
  double ne    = Getne (r);
  double cbeta = Getcbeta (r);

  double di = sqrt (M * m_p /ne /e/e /mu_0);

  return cbeta * di;
}

// ##############################
// Function to return value of QE
// ##############################
double TearX::GetQE (double r)
{
  double tauR = GettauR  (r);
  double tauA = GettauA  (r);
  double wE   = GetwE    (r);
  double s    = GetShear (r);
  double S    = tauR /tauA;

  return - pow (S, 1./3.) * ntor * wE * tauA / pow (ntor * fabs (s), 2./3.);
}

// ##############################
// Function to return value of Qe
// ##############################
double TearX::GetQe (double r)
{
  double tauR = GettauR  (r);
  double tauA = GettauA  (r);
  double we   = Getwaste (r);
  double s    = GetShear (r);
  double S    = tauR /tauA;

  return - pow (S, 1./3.) * ntor * we * tauA / pow (ntor * fabs (s), 2./3.);
}

// ##############################
// Function to return value of Qi
// ##############################
double TearX::GetQi (double r)
{
  double tauR = GettauR  (r);
  double tauA = GettauA  (r);
  double wi   = Getwasti (r);
  double s    = GetShear (r);
  double S    = tauR /tauA;

  return - pow (S, 1./3.) * ntor * wi * tauA / pow (ntor * fabs (s), 2./3.);
}

// ##############################
// Function to return value of D
// ##############################
double TearX::GetD (double r)
{
  double tauR  = GettauR  (r);
  double tauA  = GettauA  (r);
  double we    = Getwaste (r);
  double wi    = Getwasti (r);
  double dbeta = Getdbeta (r);
  double s     = GetShear (r);
  double S     = tauR /tauA;
  double iotae = - we /(wi - we);

  return pow (S * ntor * fabs (s), 1./3.) * sqrt (iotae) * dbeta /a/r;
}

// ###############################
// Function to return value of P_E
// ###############################
double TearX::GetP_E (double r)
{
  double tauR = GettauR (r);
  double tauE = GettauE (r);

  return tauR /tauE;
}

// #################################
// Function to return value of P_phi
// #################################
double TearX::GetP_phi (double r)
{
  double tauR = GettauR (r);
  double taup = Gettaup (r);

  return tauR /taup;
}

// ########################################
// Function to return value of scale factor
// ########################################
double TearX::GetScale (double r)
{
  double tauR = GettauR  (r);
  double tauA = GettauA  (r);
  double s    = GetShear (r);
  double S    = tauR /tauA;

  return pow (S * ntor * fabs (s), 1./3.);
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

