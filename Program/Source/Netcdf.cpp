// Netcdf.cpp

#include "Equilibrium.h"
#include "TJ.h"

// #################################################
// Function to write equilibrium data to netcdf file
// #################################################
void Equilibrium::WriteNetcdf (double sa, double G1, double G2)
{
  printf ("Writing data to netcdf file Plots/Equilibrium.nc:\n");

  double* Hndata  = new double[(Ns+1)*(Nr+1)];
  double* Hnpdata = new double[(Ns+1)*(Nr+1)];
  double* Vndata  = new double[(Ns+1)*(Nr+1)];
  double* Vnpdata = new double[(Ns+1)*(Nr+1)];
  double* Rdata   = new double[Nf*(Nw+1)];
  double* Zdata   = new double[Nf*(Nw+1)];
  double* rdata   = new double[Nf*(Nw+1)];
  double* tdata   = new double[Nf*(Nw+1)];
  double* Hna     = new double[Ns+1];
  double* Vna     = new double[Ns+1];
  double* npol    = new double[Ns+1];

  for (int n = 0; n <= Ns; n++)
    npol[n] = double (n);
  
  Hna[0] = 0.;
  Vna[0] = 0.;
  Hna[1] = HHfunc(1, Nr);
  Vna[1] = 0.;
  for (int n = 2; n <= Ns; n++)
    {
      Hna[n] = HHfunc(n, Nr);
      Vna[n] = VVfunc(n, Nr);
    }

  for (int n = 0; n <= Ns; n++)
    for (int i = 0; i <= Nr; i++)
      {
	Hndata [i + n*(Nr+1)] = HHfunc(n, i);
	Hnpdata[i + n*(Nr+1)] = HPfunc(n, i);
	Vndata [i + n*(Nr+1)] = VVfunc(n, i);
	Vnpdata[i + n*(Nr+1)] = VPfunc(n, i);
      }

  for (int n = 0; n < Nf; n++)
    for (int i = 0; i <= Nw; i++)
      {
	Rdata[i + n*(Nw+1)] = RR    (n, i);
	Zdata[i + n*(Nw+1)] = ZZ    (n, i);
	rdata[i + n*(Nw+1)] = rvals (n, i);
	tdata[i + n*(Nw+1)] = thvals(n, i);
      }

  try
    {
      NcFile dataFile ("Plots/Equilibrium.nc", NcFile::replace);

      NcDim p_d = dataFile.addDim ("Np", 4);
      NcDim r_d = dataFile.addDim ("Nr", Nr+1);
      NcDim s_d = dataFile.addDim ("Ns", Ns+1);
      NcDim f_d = dataFile.addDim ("Nf", Nf);
      NcDim w_d = dataFile.addDim ("Nw", Nw+1);

      double para[4];
      para[0] = epsa;
      para[1] = sa;
      para[2] = G1;
      para[3] = G2;

      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> flux_d;
      flux_d.push_back (f_d);
      flux_d.push_back (w_d);

      NcVar p_x = dataFile.addVar ("para", ncDouble, p_d);
      p_x.putVar (para);

      NcVar r_x = dataFile.addVar ("r", ncDouble, r_d);
      r_x.putVar (rr);
      NcVar g2_x = dataFile.addVar ("g_2", ncDouble, r_d);
      g2_x.putVar (g2);
      NcVar p2_x = dataFile.addVar ("p_2", ncDouble, r_d);
      p2_x.putVar (p2);
      NcVar pp_x = dataFile.addVar ("pp", ncDouble, r_d);
      pp_x.putVar (pp);
      NcVar ppp_x = dataFile.addVar ("ppp", ncDouble, r_d);
      ppp_x.putVar (ppp);
      NcVar f1_x = dataFile.addVar ("f_1", ncDouble, r_d);
      f1_x.putVar (f1);
      NcVar f3_x = dataFile.addVar ("f_3", ncDouble, r_d);
      f3_x.putVar (f3);
      NcVar q0_x = dataFile.addVar ("q_0", ncDouble, r_d);
      q0_x.putVar (q0);
      NcVar q2_x = dataFile.addVar ("q_2", ncDouble, r_d);
      q2_x.putVar (q2);
      NcVar It_x = dataFile.addVar ("I_t", ncDouble, r_d);
      It_x.putVar (It);
      NcVar Ip_x = dataFile.addVar ("I_p", ncDouble, r_d);
      Ip_x.putVar (Ip);
      NcVar Jt_x = dataFile.addVar ("J_t", ncDouble, r_d);
      Jt_x.putVar (Jt);
      NcVar Jp_x = dataFile.addVar ("J_p", ncDouble, r_d);
      Jp_x.putVar (Jp);
      NcVar q_x = dataFile.addVar ("q", ncDouble, r_d);
      q_x.putVar (q2);
      NcVar qq_x = dataFile.addVar ("qq", ncDouble, r_d);
      qq_x.putVar (qq);
      NcVar qqq_x = dataFile.addVar ("qqq", ncDouble, r_d);
      qqq_x.putVar (qqq);
      NcVar s_x = dataFile.addVar ("s", ncDouble, r_d);
      s_x.putVar (s);
      NcVar s2_x = dataFile.addVar ("s2", ncDouble, r_d);
      s2_x.putVar (s2);
      NcVar S1_x = dataFile.addVar ("S1", ncDouble, r_d);
      S1_x.putVar (S1);
      NcVar S2_x = dataFile.addVar ("S2", ncDouble, r_d);
      S2_x.putVar (S2);
      NcVar P1_x = dataFile.addVar ("P1", ncDouble, r_d);
      P1_x.putVar (P1);
      NcVar P2_x = dataFile.addVar ("P2", ncDouble, r_d);
      P2_x.putVar (P2);
      NcVar P3_x = dataFile.addVar ("P3", ncDouble, r_d);
      P3_x.putVar (P3);
      NcVar P3a_x = dataFile.addVar ("P3a", ncDouble, r_d);
      P3a_x.putVar (P3a);
 
      NcVar Hn_x = dataFile.addVar ("Hn", ncDouble, shape_d);
      Hn_x.putVar (Hndata);
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (Hnpdata);
      NcVar Vn_x = dataFile.addVar ("Vn", ncDouble, shape_d);
      Vn_x.putVar (Vndata);
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (Vnpdata);

      NcVar R_x = dataFile.addVar ("R", ncDouble, flux_d);
      R_x.putVar (Rdata);
      NcVar Z_x = dataFile.addVar ("Z", ncDouble, flux_d);
      Z_x.putVar (Zdata);
      NcVar rr_x = dataFile.addVar ("rr", ncDouble, flux_d);
      rr_x.putVar (rdata);
      NcVar t_x = dataFile.addVar ("theta", ncDouble, flux_d);
      t_x.putVar (tdata);

      NcVar n_x = dataFile.addVar ("n", ncDouble, s_d);
      n_x.putVar (npol);
      NcVar Hna_x = dataFile.addVar ("Hna", ncDouble, s_d);
      Hna_x.putVar (Hna);
      NcVar Vna_x = dataFile.addVar ("Vna", ncDouble, s_d);
      Vna_x.putVar (Vna);

      NcVar Rbound_x = dataFile.addVar ("Rbound", ncDouble, w_d);
      Rbound_x.putVar (Rbound);
      NcVar Zbound_x = dataFile.addVar ("Zbound", ncDouble, w_d);
      Zbound_x.putVar (Zbound);
      NcVar tbound_x = dataFile.addVar ("tbound", ncDouble, w_d);
      tbound_x.putVar (tbound);
      NcVar wbound_x = dataFile.addVar ("wbound", ncDouble, w_d);
      wbound_x.putVar (wbound);
      NcVar R2b_x = dataFile.addVar ("R2bound", ncDouble, w_d);
      R2b_x.putVar (R2b);
      NcVar grr2b_x = dataFile.addVar ("grr2bound", ncDouble, w_d);
      grr2b_x.putVar (grr2b);
      NcVar dRdtheta_x = dataFile.addVar ("dRdtheta", ncDouble, w_d);
      dRdtheta_x.putVar (dRdtheta);
      NcVar dZdtheta_x = dataFile.addVar ("dZdtheta", ncDouble, w_d);
      dZdtheta_x.putVar (dZdtheta);
    }
  catch(NcException& e)
    {
      printf ("%s\n", e.what ());
      printf ("Error writing data to netcdf file Plots/Equilbrium.nc\n");
      exit (1);
    }

  delete[] Hndata; delete[] Hnpdata; delete[] Vndata; delete[] Vnpdata;
  delete[] Rdata;  delete[] Zdata;   delete[] Hna;    delete[] Vna;
  delete[] rdata;  delete[] tdata;   delete[] npol;
}

// ##################################################
// Function to read equilibrium data from netcdf file
// ##################################################
void TJ::ReadNetcdf ()
{
  printf ("Reading data from netcdf file Plots/Equilibrium.nc:\n");
  
  try
    {
      NcFile dataFile ("Plots/Equilibrium.nc", NcFile::read);
      
      NcVar p_x = dataFile.getVar ("para");
      NcDim p_d = p_x.getDim (0);

      int     Np   = p_d.getSize ();
      double* para = new double[Np];

      p_x.getVar(para);
      epsa = para[0];
      sa   = para[1];
      G1   = para[2];
      G2   = para[3];
      
      NcVar r_x   = dataFile.getVar ("r");
      NcVar pp_x  = dataFile.getVar ("pp");
      NcVar ppp_x = dataFile.getVar ("ppp");
      NcVar q_x   = dataFile.getVar ("q");
      NcVar s_x   = dataFile.getVar ("s");
      NcVar s2_x  = dataFile.getVar ("s2");
      NcVar S1_x  = dataFile.getVar ("S1");
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P3_x  = dataFile.getVar ("P3");
      NcDim r_d   = r_x.getDim (0);

      Nr  = r_d.getSize () - 1;
      rr  = new double[Nr+1];
      pp  = new double[Nr+1];
      ppp = new double[Nr+1];
      q   = new double[Nr+1];
      s   = new double[Nr+1];
      s2  = new double[Nr+1];
      S1  = new double[Nr+1];
      P1  = new double[Nr+1];
      P2  = new double[Nr+1];
      P3  = new double[Nr+1];

      r_x.  getVar (rr);
      pp_x. getVar (pp);
      ppp_x.getVar (ppp);
      q_x  .getVar (q);
      s_x  .getVar (s);
      s2_x .getVar (s2);
      S1_x .getVar (S1);
      P1_x .getVar (P1);
      P2_x .getVar (P2);
      P3_x .getVar (P3);

      NcVar Hn_x  = dataFile.getVar ("Hn");
      NcVar Hnp_x = dataFile.getVar ("Hnp");
      NcVar Vn_x  = dataFile.getVar ("Vn");
      NcVar Vnp_x = dataFile.getVar ("Vnp");
      NcDim s_d   = Hn_x.getDim (0);

      Ns              = s_d.getSize () - 1;
      double* Hndata  = new double[(Ns+1)*(Nr+1)];
      double* Hnpdata = new double[(Ns+1)*(Nr+1)];
      double* Vndata  = new double[(Ns+1)*(Nr+1)];
      double* Vnpdata = new double[(Ns+1)*(Nr+1)];

      Hn_x .getVar (Hndata);
      Hnp_x.getVar (Hnpdata);
      Vn_x .getVar (Vndata);
      Vnp_x.getVar (Vnpdata);

      HHfunc.resize (Ns+1, Nr+1);
      VVfunc.resize (Ns+1, Nr+1);
      HPfunc.resize (Ns+1, Nr+1);
      VPfunc.resize (Ns+1, Nr+1);
      
      for (int n = 0; n <= Ns; n++)
	for (int i = 0; i <= Nr; i++)
	  {
	    HHfunc(n, i) = Hndata [i + n*(Nr+1)];
	    HPfunc(n, i) = Hnpdata[i + n*(Nr+1)];
	    VVfunc(n, i) = Vndata [i + n*(Nr+1)];
	    VPfunc(n, i) = Vnpdata[i + n*(Nr+1)];
	  }

      NcVar RR_x = dataFile.getVar ("R");
      NcVar ZZ_x = dataFile.getVar ("Z");
      NcVar rr_x = dataFile.getVar ("rr");
      NcVar tt_x = dataFile.getVar ("theta");
      NcDim f_d  = RR_x.getDim (0);
      NcDim w_d  = RR_x.getDim (1);

      Nf = f_d.getSize ();
      Nw = w_d.getSize () - 1;

      RR.resize     (Nf, Nw+1);
      ZZ.resize     (Nf, Nw+1);
      rvals.resize  (Nf, Nw+1);
      thvals.resize (Nf, Nw+1);

      double* RRdata = new double[Nf*(Nw+1)];
      double* ZZdata = new double[Nf*(Nw+1)];
      double* rrdata = new double[Nf*(Nw+1)];
      double* ttdata = new double[Nf*(Nw+1)];

      RR_x.getVar (RRdata);
      ZZ_x.getVar (ZZdata);
      rr_x.getVar (rrdata);
      tt_x.getVar (ttdata);

      rf = new double[Nf];
      for (int n = 0; n < Nf; n++)
	rf[n] = rrdata[n*Nw];
      
      for (int n = 0; n < Nf; n++)
	for (int i = 0; i <= Nw; i++)
	  {
	    RR    (n, i) = RRdata[i + n*(Nw+1)];
	    ZZ    (n, i) = ZZdata[i + n*(Nw+1)];
	    rvals (n, i) = rrdata[i + n*(Nw+1)];
	    thvals(n, i) = ttdata[i + n*(Nw+1)];
	  }

      NcVar tbound_x = dataFile.getVar ("tbound");
      NcVar Rbound_x = dataFile.getVar ("Rbound");
      NcVar Zbound_x = dataFile.getVar ("Zbound");
      NcVar dRdthe_x = dataFile.getVar ("dRdtheta");
      NcVar dZdthe_x = dataFile.getVar ("dZdtheta");
 
      tbound = new double[Nw+1];
      Rbound = new double[Nw+1];
      Zbound = new double[Nw+1];
      dRdthe = new double[Nw+1];
      dZdthe = new double[Nw+1];

      tbound_x.getVar (tbound);
      Rbound_x.getVar (Rbound);
      Zbound_x.getVar (Zbound);
      dRdthe_x.getVar (dRdthe);
      dZdthe_x.getVar (dZdthe);
      
      delete[] para;   delete[] Hndata; delete[] Hnpdata; delete[] Vndata; delete[] Vnpdata;
      delete[] RRdata; delete[] ZZdata; delete[] rrdata;  delete[] ttdata;
    }
  catch (NcException& e)
     {
       printf ("%s\n", e.what ());
       printf ("Error reading data from netcdf file Inputs/Equilbrium.nc\n");
       exit (1);
     }
}

// ###############################################
// Function to write stability data to netcdf file
// ###############################################
void TJ::WriteNetcdf ()
{
  printf ("Writing data to netcdf file TJ.cpp:\n");
 
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
      printf ("%s\n", e.what ());
      printf ("Error writing data to netcdf file Plots/TJ.nc\n");
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