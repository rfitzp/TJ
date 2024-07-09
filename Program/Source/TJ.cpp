// TJ.cpp

#include "TJ.h"

// ###########
// Constructor
// ###########
TJ::TJ ()
{
  // ---------------------------
  // Set root finding parameters
  // ---------------------------
  Eta     = 1.e-16;
  Maxiter = 30;
  
  // -----------------------------------
  // Set adaptive integration parameters
  // -----------------------------------
  maxrept = 50;
  flag    = 2;
  
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

  // --------------------------------
  // Read namelist file Inputs/TJ.nml
  // --------------------------------
  NameListTJ (&NTOR, &MMIN, &MMAX, 
	      &EPS, &DEL, &NFIX, &NDIAG, &NULC, &ITERMAX,
	      &FREE, &SYMM, 
	      &acc, &h0, &hmin, &hmax, &EPSF);

  // ............
  // Sanity check
  // ............
  if (NTOR < 1)
    {
      printf ("TJ: Error - NTOR must be positive");
      exit (1);
    }
  if (MMAX < MMIN)
    {
      printf ("TJ: Error - MMIN must be less that MMAX");
      exit (1);
    }
  if (MMAX < MMIN)
    {
      printf ("TJ: Error - MMIN must be less that MMAX");
      exit (1);
    }
  if (EPS <= 0.)
    {
      printf ("TJ: Error - EPS must be positive");
      exit (1);
    }
  if (DEL <= 0.)
    {
      printf ("TJ: Error - DEL must be positive");
      exit (1);
    }
  if (NFIX < 0)
    {
      printf ("TJ: Error - NFIX cannot be negative");
      exit (1);
    }
  if (NDIAG < 0)
    {
      printf ("TJ: Error - NDIAG cannot be less that two");
      exit (1);
    }
  if (NULC <= 0.)
    {
      printf ("TJ: Error - NULC must be positive");
      exit (1);
    }
  if (ITERMAX < 0)
    {
      printf ("TJ: Error - ITERMAX cannot be negative");
      exit (1);
    }
    if (acc <= 0.)
    {
      printf ("TJ:: Error - acc must be positive\n");
      exit (1);
    }
  if (h0 <= 0.)
    {
      printf ("TJ:: Error - h0 must be positive\n");
      exit (1);
    }
  if (hmin <= 0.)
    {
      printf ("TJ:: Error - hmin must be positive\n");
      exit (1);
    }
  if (hmax <= 0.)
    {
      printf ("TJ:: Error - hmax must be positive\n");
      exit (1);
    }
  if (hmax < hmin)
    {
      printf ("TJ:: Error - hmax must exceed hmin\n");
      exit (1);
    }
  if (EPSF <= 0.)
    {
      printf ("TJ: Error - EPSF must be positive");
      exit (1);
    }
  
  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\nSubprogram TJ::\n");
  printf ("Calculation parameters:\n");
  printf ("ntor = %3d        mmin  = %3d        mmax = %3d        eps     = %10.3e del  = %10.3e\n",
	  NTOR, MMIN, MMAX, EPS, DEL);
  printf ("nfix = %3d        ndiag = %3d       nulc = %10.3e itermax = %3d        free =  %1d         symm = %1d\n",
	  NFIX, NDIAG, NULC, ITERMAX, FREE, SYMM);
  printf ("acc  = %10.3e h0    = %10.3e hmin = %10.3e hmax    = %10.3e epsf = %10.3e\n",
	  acc, h0, hmin, hmax, EPSF);
}

// ##########
// Destructor
// ##########
TJ::~TJ ()
{
}

// #########################
// Function to solve problem
// #########################
void TJ::Solve ()
{ 
  // Set toroidal and poloidal mode numbers
  SetModeNumbers ();

  // Read equlibrium data
  ReadEquilibrium ();

  // Calculate metric data at plasma boundary
  CalculateMetric ();

  // Calculate vacuum matrices
  GetVacuum ();

  // Find rational surfaces
  FindRational ();

  // Solve outer region odes
  ODESolve ();

  // Determine tearing mode dispersion relation and tearing eigenfunctions
  FindDispersion ();
  
  // Write program data
  printf ("Writing data to netcdf file:\n");
  WriteNetCDF ();

  // Clean up
  CleanUp ();
}

// ##################################################
// Function to set toroidal and poloidal mode numbers
// ##################################################
void TJ::SetModeNumbers ()
{
  ntor = double (NTOR);
  J    = MMAX - MMIN + 1;
  MPOL = new int[J];
  mpol = new double[J];

  for (int j = 0; j < J; j++)
    {
      MPOL[j] = MMIN + j;
      mpol[j] = double (MMIN + j);
    }
}

// #######################################
// Function to write data to Outputs/TJ.nc
// #######################################
void TJ::WriteNetCDF ()
{
  // ..........................
  // Output data to netcdf file
  // ..........................
  double* Hndata  = new double[(Ns+1)*(Nr+1)];
  double* Hnpdata = new double[(Ns+1)*(Nr+1)];
  double* Vndata  = new double[(Ns+1)*(Nr+1)];
  double* Vnpdata = new double[(Ns+1)*(Nr+1)];
  double* Lmmp_r  = new double[(Nr+1)*J*J];
  double* Mmmp_r  = new double[(Nr+1)*J*J];
  double* Nmmp_r  = new double[(Nr+1)*J*J];
  double* Pmmp_r  = new double[(Nr+1)*J*J];
  double* Lmmp_i  = new double[(Nr+1)*J*J];
  double* Mmmp_i  = new double[(Nr+1)*J*J];
  double* Nmmp_i  = new double[(Nr+1)*J*J];
  double* Pmmp_i  = new double[(Nr+1)*J*J];
  double* Ltest   = new double[(Nr+1)*J*J];
  double* MNtest  = new double[(Nr+1)*J*J];
  double* Ptest   = new double[(Nr+1)*J*J];
  double* Pvac_r  = new double[J*J];
  double* Pvac_i  = new double[J*J];
  double* Rvac_r  = new double[J*J];
  double* Rvac_i  = new double[J*J];
  double* Avac_r  = new double[J*J];
  double* Avac_i  = new double[J*J];
  double* Hmat_r  = new double[J*J];
  double* Hmat_i  = new double[J*J];
  double* Hres_r  = new double[J*J];
  double* Hres_i  = new double[J*J];
  double* Ttest_i = new double[K*NDIAG];
  double* Pnorm_i = new double[K*NDIAG];
  double* Znorm_i = new double[K*NDIAG];
  double* PPPsi_r = new double[J*K*NDIAG];
  double* PPPsi_i = new double[J*K*NDIAG];
  double* ZZZ_r   = new double[J*K*NDIAG];
  double* ZZZ_i   = new double[J*K*NDIAG];
  double* PPF_r   = new double[J*nres*NDIAG];
  double* PPF_i   = new double[J*nres*NDIAG];
  double* ZZF_r   = new double[J*nres*NDIAG];
  double* ZZF_i   = new double[J*nres*NDIAG];
  double* PPU_r   = new double[J*nres*NDIAG];
  double* PPU_i   = new double[J*nres*NDIAG];
  double* ZZU_r   = new double[J*nres*NDIAG];
  double* ZZU_i   = new double[J*nres*NDIAG];
  double* TTf     = new double[nres*NDIAG];
  double* TTu     = new double[nres*NDIAG];
  double* TFull   = new double[nres*nres*NDIAG];
  double* TUnrc   = new double[nres*nres*NDIAG];
  double* PPV_r   = new double[nres*Nf*(Nw+1)];
  double* PPV_i   = new double[nres*Nf*(Nw+1)];
  double* ZZV_r   = new double[nres*Nf*(Nw+1)];
  double* ZZV_i   = new double[nres*Nf*(Nw+1)];

  for (int n = 0; n <= Ns; n++)
    for (int i = 0; i <= Nr; i++)
      {
	Hndata [i + n*(Nr+1)] = HHfunc(n, i);
	Hnpdata[i + n*(Nr+1)] = HPfunc(n, i);
	Vndata [i + n*(Nr+1)] = VVfunc(n, i);
	Vnpdata[i + n*(Nr+1)] = VPfunc(n, i);
      }

  int cnt = 0;
  for (int i = 0; i <= Nr; i++)
    for (int j = 0; j < J; j++)
      for (int jp = 0; jp < J; jp++)
	{
	  complex<double> Lmmp, Mmmp, Nmmp, Pmmp;
	  complex<double> Lmpm, Mmpm, Nmpm, Pmpm;
	  GetMatrices (rr[i], MPOL[j], MPOL[jp], Lmmp, Mmmp, Nmmp, Pmmp);
	  GetMatrices (rr[i], MPOL[jp], MPOL[j], Lmpm, Mmpm, Nmpm, Pmpm);

	  Lmmp_r[cnt] = Lmmp.real();
	  Lmmp_i[cnt] = Lmmp.imag();
	  Mmmp_r[cnt] = Mmmp.real();
	  Mmmp_i[cnt] = Mmmp.imag();
	  Nmmp_r[cnt] = Nmmp.real();
	  Nmmp_i[cnt] = Nmmp.imag();
	  Pmmp_r[cnt] = Pmmp.real();
	  Pmmp_i[cnt] = Pmmp.imag();

	  Ltest[cnt]  = abs (Lmmp - conj (Lmpm));
	  MNtest[cnt] = abs (Mmmp + conj (Nmpm));
	  Ptest[cnt]  = abs (Pmmp - conj (Pmpm));
	  cnt++;
	}

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Pvac_r[cnt] = Pvac(j, jp).real();
	Pvac_i[cnt] = Pvac(j, jp).imag();
	Rvac_r[cnt] = Rvac(j, jp).real();
	Rvac_i[cnt] = Rvac(j, jp).imag();
	Avac_r[cnt] = Avac(j, jp).real();
	Avac_i[cnt] = Avac(j, jp).imag();
	Hmat_r[cnt] = Hmat(j, jp).real();
	Hmat_i[cnt] = Hmat(j, jp).imag();
	Hres_r[cnt] = Hmat(j, jp).real() - Hdag(j, jp).real();
	Hres_i[cnt] = Hmat(j, jp).imag() - Hdag(j, jp).imag();
	cnt++;
      }

  cnt = 0;
  for (int j = 0; j < K; j++)
    for (int i = 0; i < NDIAG; i++)
      {
	Ttest_i[cnt] = Ttest (j, i);
	Pnorm_i[cnt] = Pnorm (j, i);
	Znorm_i[cnt] = Znorm (j, i);
	cnt++;
      }

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < K; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  PPPsi_r[cnt] = real (YYY(j,   jp, i));
	  PPPsi_i[cnt] = imag (YYY(j,   jp, i));
	  ZZZ_r  [cnt] = real (YYY(J+j, jp, i));
	  ZZZ_i  [cnt] = imag (YYY(J+j, jp, i));
	  cnt++;
	}

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < nres; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  PPF_r[cnt] = real (Psif(j, jp, i));
	  PPF_i[cnt] = imag (Psif(j, jp, i));
	  ZZF_r[cnt] = real (Zf  (j, jp, i));
	  ZZF_i[cnt] = imag (Zf  (j, jp, i));
	  cnt++;
	}

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < nres; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  PPU_r[cnt] = real (Psiu(j, jp, i));
	  PPU_i[cnt] = imag (Psiu(j, jp, i));
	  ZZU_r[cnt] = real (Zu  (j, jp, i));
	  ZZU_i[cnt] = imag (Zu  (j, jp, i));
	  cnt++;
	}

  cnt = 0;
  for (int jp = 0; jp < nres; jp++)
    for (int i = 0; i < NDIAG; i++)
      {
	TTf[cnt] = Tf(jp, i);
	TTu[cnt] = Tu(jp, i);
	cnt++;
      }

  cnt = 0;
  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  TFull[cnt] = Tfull(j, jp, i);
	  TUnrc[cnt] = Tunrc(j, jp, i);
	  cnt++;
	}

  cnt = 0;
  for (int k = 0; k < nres; k++)
    for (int i = 0; i < Nf; i++)
      for (int l = 0; l <= Nw; l++)
	{
	  PPV_r[cnt] = real (Psiuv(k, i, l));
	  PPV_i[cnt] = imag (Psiuv(k, i, l));
	  ZZV_r[cnt] = real (Zuv  (k, i, l));
	  ZZV_i[cnt] = imag (Zuv  (k, i, l));
	  cnt++;
	}

  try
    {
      NcFile dataFile ("Plots/TJ.nc", NcFile::replace);

      NcDim r_d = dataFile.addDim ("Nr",    Nr+1);
      NcDim s_d = dataFile.addDim ("Ns",    Ns+1);
      NcDim x_d = dataFile.addDim ("nres",  nres);
      NcDim j_d = dataFile.addDim ("J",     J);
      NcDim k_d = dataFile.addDim ("K",     K);
      NcDim d_d = dataFile.addDim ("ndiag", NDIAG);
      NcDim f_d = dataFile.addDim ("Nf",    Nf);
      NcDim w_d = dataFile.addDim ("Nw",    Nw+1);

      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> matrix_d;
      matrix_d.push_back (r_d);
      matrix_d.push_back (j_d);
      matrix_d.push_back (j_d);

      vector<NcDim> vacuum_d;
      vacuum_d.push_back (j_d);
      vacuum_d.push_back (j_d);

      vector<NcDim> torque_d;
      torque_d.push_back (k_d);
      torque_d.push_back (d_d);

      vector<NcDim> chi_d;
      chi_d.push_back (x_d);
      chi_d.push_back (j_d);
      
      vector<NcDim> soln_d;
      soln_d.push_back (j_d);
      soln_d.push_back (k_d);
      soln_d.push_back (d_d);

      vector<NcDim> full_d;
      full_d.push_back (j_d);
      full_d.push_back (x_d);
      full_d.push_back (d_d);

      vector<NcDim> t_d;
      t_d.push_back (x_d);
      t_d.push_back (d_d);

      vector<NcDim> tt_d;
      tt_d.push_back (x_d);
      tt_d.push_back (x_d);
      tt_d.push_back (d_d);

      vector<NcDim> v_d;
      v_d.push_back (x_d);
      v_d.push_back (f_d);
      v_d.push_back (w_d);
      
      NcVar r_x = dataFile.addVar ("r", ncDouble, r_d);
      r_x.putVar (rr);
      NcVar pp_x = dataFile.addVar ("pp", ncDouble, r_d);
      pp_x.putVar (pp);
      NcVar ppp_x = dataFile.addVar ("ppp", ncDouble, r_d);
      ppp_x.putVar (ppp);
      NcVar q_x = dataFile.addVar ("q", ncDouble, r_d);
      q_x.putVar (q);
      NcVar s_x = dataFile.addVar ("s", ncDouble, r_d);
      s_x.putVar (s);
      NcVar s2_x = dataFile.addVar ("s2", ncDouble, r_d);
      s2_x.putVar (s2);
      NcVar S1_x = dataFile.addVar ("S1", ncDouble, r_d);
      S1_x.putVar (S1);
      NcVar P1_x = dataFile.addVar ("P1", ncDouble, r_d);
      P1_x.putVar (P1);
      NcVar P2_x = dataFile.addVar ("P2", ncDouble, r_d);
      P2_x.putVar (P2);
      NcVar P3_x = dataFile.addVar ("P3", ncDouble, r_d);
      P3_x.putVar (P3);
 
      NcVar Hn_x = dataFile.addVar ("Hn", ncDouble, shape_d);
      Hn_x.putVar (Hndata);
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (Hnpdata);
      NcVar Vn_x = dataFile.addVar ("Vn", ncDouble, shape_d);
      Vn_x.putVar (Vndata);
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (Vnpdata);

      NcVar rres_x = dataFile.addVar ("rres", ncDouble, x_d);
      rres_x.putVar (rres);

      NcVar mpol_x = dataFile.addVar ("mpol", ncInt, j_d);
      mpol_x.putVar (MPOL);

      NcVar lmmpr_x = dataFile.addVar ("Lmmp_r", ncDouble, matrix_d);
      lmmpr_x.putVar (Lmmp_r);
      NcVar lmmpi_x = dataFile.addVar ("Lmmp_i", ncDouble, matrix_d);
      lmmpi_x.putVar (Lmmp_i);

      NcVar mmmpr_x = dataFile.addVar ("Mmmp_r", ncDouble, matrix_d);
      mmmpr_x.putVar (Mmmp_r);
      NcVar mmmpi_x = dataFile.addVar ("Mmmp_i", ncDouble, matrix_d);
      mmmpi_x.putVar (Mmmp_i);

      NcVar nmmpr_x = dataFile.addVar ("Nmmp_r", ncDouble, matrix_d);
      nmmpr_x.putVar (Nmmp_r);
      NcVar nmmpi_x = dataFile.addVar ("Nmmp_i", ncDouble, matrix_d);
      nmmpi_x.putVar (Nmmp_i);

      NcVar pmmpr_x = dataFile.addVar ("Pmmp_r", ncDouble, matrix_d);
      pmmpr_x.putVar (Pmmp_r);
      NcVar pmmpi_x = dataFile.addVar ("Pmmp_i", ncDouble, matrix_d);
      pmmpi_x.putVar (Pmmp_i);

      NcVar ltest_x = dataFile.addVar ("Ltest", ncDouble, matrix_d);
      ltest_x.putVar (Ltest);
      NcVar mntest_x = dataFile.addVar ("MNtest", ncDouble, matrix_d);
      mntest_x.putVar (MNtest);
      NcVar ptest_x = dataFile.addVar ("Ptest", ncDouble, matrix_d);
      ptest_x.putVar (Ptest);

      NcVar pvacr_x = dataFile.addVar ("Pvac_r", ncDouble, vacuum_d);
      pvacr_x.putVar (Pvac_r);
      NcVar pvaci_x = dataFile.addVar ("Pvac_i", ncDouble, vacuum_d);
      pvaci_x.putVar (Pvac_i);
      NcVar rvacr_x = dataFile.addVar ("Rvac_r", ncDouble, vacuum_d);
      rvacr_x.putVar (Rvac_r);
      NcVar rvaci_x = dataFile.addVar ("Rvac_i", ncDouble, vacuum_d);
      rvaci_x.putVar (Rvac_i);
      
      NcVar avacr_x = dataFile.addVar ("Avac_r", ncDouble, vacuum_d);
      avacr_x.putVar (Avac_r);
      NcVar avaci_x = dataFile.addVar ("Avac_i", ncDouble, vacuum_d);
      avaci_x.putVar (Avac_i);
 
      NcVar hmatr_x = dataFile.addVar ("Hmat_r", ncDouble, vacuum_d);
      hmatr_x.putVar (Hmat_r);
      NcVar hmati_x = dataFile.addVar ("Hmat_i", ncDouble, vacuum_d);
      hmati_x.putVar (Hmat_i);
      NcVar hresr_x = dataFile.addVar ("Hres_r", ncDouble, vacuum_d);
      hresr_x.putVar (Hres_r);
      NcVar hresi_x = dataFile.addVar ("Hres_i", ncDouble, vacuum_d);
      hresi_x.putVar (Hres_i);
 
      NcVar rgrid_x = dataFile.addVar ("r_grid", ncDouble, d_d);
      rgrid_x.putVar (Rgrid);
      NcVar hode_x = dataFile.addVar ("h_ode", ncDouble, d_d);
      hode_x.putVar (hode);
      NcVar eode_x = dataFile.addVar ("err_ode", ncDouble, d_d);
      eode_x.putVar (eode);

      NcVar t_x = dataFile.addVar ("theta", ncDouble, w_d);
      t_x.putVar (tbound);
      NcVar cmu_x = dataFile.addVar ("cosmu", ncDouble, w_d);
      cmu_x.putVar (cmu);
      NcVar e_x = dataFile.addVar ("eta", ncDouble, w_d);
      e_x.putVar (eeta);
      NcVar ceta_x = dataFile.addVar ("coseta", ncDouble, w_d);
      ceta_x.putVar (ceta);
      NcVar seta_x = dataFile.addVar ("sineta", ncDouble, w_d);
      seta_x.putVar (seta);
      NcVar R2grgz_x = dataFile.addVar ("R2grgz", ncDouble, w_d);
      R2grgz_x.putVar (R2grgz);
      NcVar R2grge_x = dataFile.addVar ("R2grge", ncDouble, w_d);
      R2grge_x.putVar (R2grge);

      NcVar ttest_x = dataFile.addVar ("Torque_test", ncDouble, torque_d);
      ttest_x.putVar (Ttest_i);
      NcVar pnorm_x = dataFile.addVar ("Psi_norm", ncDouble, torque_d);
      pnorm_x.putVar (Pnorm_i);
      NcVar znorm_x = dataFile.addVar ("Z_norm", ncDouble, torque_d);
      znorm_x.putVar (Znorm_i);

      NcVar pppsir_x = dataFile.addVar ("Psi_r", ncDouble, soln_d);
      pppsir_x.putVar (PPPsi_r);
      NcVar pppsii_x = dataFile.addVar ("Psi_i", ncDouble, soln_d);
      pppsii_x.putVar (PPPsi_i);
      NcVar zzzr_x   = dataFile.addVar ("Z_r", ncDouble, soln_d);
      zzzr_x.putVar (ZZZ_r);
      NcVar zzzi_x   = dataFile.addVar ("Z_i", ncDouble, soln_d);
      zzzi_x.putVar (ZZZ_i);

      NcVar ppfr_x = dataFile.addVar ("Psi_full_r", ncDouble, full_d);
      ppfr_x.putVar (PPF_r);
      NcVar ppfi_x = dataFile.addVar ("Psi_full_i", ncDouble, full_d);
      ppfi_x.putVar (PPF_i);
      NcVar zzfr_x = dataFile.addVar ("Z_full_r", ncDouble, full_d);
      zzfr_x.putVar (ZZF_r);
      NcVar zzfi_x = dataFile.addVar ("Z_full_i", ncDouble, full_d);
      zzfi_x.putVar (ZZF_i);

      NcVar ppur_x = dataFile.addVar ("Psi_unrc_r", ncDouble, full_d);
      ppur_x.putVar (PPU_r);
      NcVar ppui_x = dataFile.addVar ("Psi_unrc_i", ncDouble, full_d);
      ppui_x.putVar (PPU_i);
      NcVar zzur_x = dataFile.addVar ("Z_unrc_r", ncDouble, full_d);
      zzur_x.putVar (ZZU_r);
      NcVar zzui_x = dataFile.addVar ("Z_unrc_i", ncDouble, full_d);
      zzui_x.putVar (ZZU_i);

      NcVar tf_x = dataFile.addVar ("Torque_full", ncDouble, t_d);
      tf_x.putVar (TTf);
      NcVar tu_x = dataFile.addVar ("Torque_unrc", ncDouble, t_d);
      tu_x.putVar (TTu);

      NcVar mres_x = dataFile.addVar ("m_res", ncInt, x_d);
      mres_x.putVar (mres);

      NcVar tfull_x = dataFile.addVar ("Torque_pair_full", ncDouble, tt_d);
      tfull_x.putVar (TFull);
      NcVar tunrc_x = dataFile.addVar ("Torque_pair_unrc", ncDouble, tt_d);
      tunrc_x.putVar (TUnrc);

      NcVar ppvr_x = dataFile.addVar ("Psi_unrc_eig_r", ncDouble, v_d);
      ppvr_x.putVar (PPV_r);
      NcVar ppvi_x = dataFile.addVar ("Psi_unrc_eig_i", ncDouble, v_d);
      ppvi_x.putVar (PPV_i);
      NcVar zzvr_x = dataFile.addVar ("Z_unrc_eig_r", ncDouble, v_d);
      zzvr_x.putVar (ZZV_r);
      NcVar zzvi_x = dataFile.addVar ("Z_unrc_eig_i", ncDouble, v_d);
      zzvi_x.putVar (ZZV_i);
    }
  catch (NcException& e)
    {
      e.what ();
      printf ("Error writing Plots/TJ.nc\n");
      exit (1);
    }

  delete[] Hndata;  delete[] Hnpdata; delete[] Vndata;  delete[] Vnpdata;
  delete[] Lmmp_r;  delete[] Mmmp_r;  delete[] Nmmp_r;  delete[] Pmmp_r;
  delete[] Lmmp_i;  delete[] Mmmp_i;  delete[] Nmmp_i;  delete[] Pmmp_i;
  delete[] Ltest;   delete[] MNtest;  delete[] Ptest;
  delete[] Avac_r;  delete[] Avac_i;  delete[] Pvac_r;  delete[] Pvac_i;
  delete[] Rvac_r;  delete[] Rvac_i;  delete[] Hmat_r;  delete[] Hmat_i;
  delete[] Hres_r;  delete[] Hres_i;  
  delete[] Ttest_i; delete[] Pnorm_i; delete[] Znorm_i;
  delete[] PPPsi_r; delete[] PPPsi_i; delete[] ZZZ_r;   delete[] ZZZ_i;
  delete[] PPF_r;   delete[] PPF_i;   delete[] ZZF_r;   delete[] ZZF_i;
  delete[] TTf;     delete[] TTu;     delete[] TFull;   delete[] TUnrc;
  delete[] PPV_r;   delete[] PPV_i;   delete[] ZZV_r;   delete[] ZZV_i;
  delete[] R2grgz;  delete[] R2grge;  delete[] cmu;     delete[] ceta;
  delete[] seta;    delete[] eeta;
}

// #############################
// Function to deallocate memory
// #############################
void TJ::CleanUp ()
{
  delete[] MPOL; delete[] mpol;
  
  delete[] rr; delete[] pp; delete[] ppp; delete[] q; delete[] s; delete[] s2;
  delete[] S1; delete[] P1; delete[] P2;  delete[] P3;

  gsl_spline_free (ppspline);
  gsl_spline_free (pppspline);
  gsl_spline_free (qspline);
  gsl_spline_free (sspline);
  gsl_spline_free (s2spline);
  gsl_spline_free (S1spline);
  gsl_spline_free (P1spline);
  gsl_spline_free (P2spline);
  gsl_spline_free (P3spline);

  gsl_interp_accel_free (ppacc);
  gsl_interp_accel_free (pppacc);
  gsl_interp_accel_free (qacc);
  gsl_interp_accel_free (sacc);
  gsl_interp_accel_free (s2acc);
  gsl_interp_accel_free (S1acc);
  gsl_interp_accel_free (P1acc);
  gsl_interp_accel_free (P2acc);
  gsl_interp_accel_free (P3acc);

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

  delete[] mres;   delete[] qres;   delete[] rres; delete[] sres; delete[] DIres;
  delete[] nuLres; delete[] nuSres; delete[] qerr; delete[] Jres;

  delete[] hode; delete[] eode; delete[] Rgrid; delete[] rf;

  gsl_spline_free (Rrzspline); gsl_spline_free (Rrespline); gsl_interp_accel_free (Rrzacc); gsl_interp_accel_free (Rreacc);
  gsl_spline_free (Rbspline);  gsl_spline_free (Zbspline);  gsl_interp_accel_free (Rbacc);  gsl_interp_accel_free (Zbacc);

  delete[] Rbound; delete[] Zbound; delete[] tbound; delete[] dRdthe; delete[] dZdthe;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* TJ::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to open file for reading
// #################################
FILE* TJ::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

