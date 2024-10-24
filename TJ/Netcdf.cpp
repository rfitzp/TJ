// Netcdf.cpp

#include "TJ.h"

// ##################################################
// Function to read equilibrium data from netcdf file
// ##################################################
void TJ::ReadNetcdf ()
{
  printf ("Reading data from netcdf file Outputs/Equilibrium/Equilibrium.nc:\n");
  
  try
    {
      NcFile dataFile ("../Outputs/Equilibrium/Equilibrium.nc", NcFile::read);
      
      NcVar p_x = dataFile.getVar ("para");
      NcDim p_d = p_x.getDim (0);

      int     Np   = p_d.getSize ();
      double* para = new double[Np];

      p_x.getVar(para);
      epsa = para[0];
      sa   = para[1];
      
      NcVar r_x   = dataFile.getVar ("r");
      NcVar P_x   = dataFile.getVar ("PsiN");
      NcVar PP_x  = dataFile.getVar ("Psi");
      NcVar f_x   = dataFile.getVar ("f");
      NcVar g2_x  = dataFile.getVar ("g_2");
      NcVar p2_x  = dataFile.getVar ("p_2");
      NcVar pp_x  = dataFile.getVar ("pp");
      NcVar ppp_x = dataFile.getVar ("ppp");
      NcVar q_x   = dataFile.getVar ("q");
      NcVar s_x   = dataFile.getVar ("s");
      NcVar s2_x  = dataFile.getVar ("s2");
      NcVar s0_x  = dataFile.getVar ("s0");
      NcVar S1_x  = dataFile.getVar ("S1");
      NcVar S3_x  = dataFile.getVar ("S3");
      NcVar S4_x  = dataFile.getVar ("S4");
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P1a_x = dataFile.getVar ("P1a");
      NcVar P2a_x = dataFile.getVar ("P2a");
      NcVar P3_x  = dataFile.getVar ("P3");
      NcVar P4_x  = dataFile.getVar ("P4");
      NcDim r_d   = r_x.getDim (0);

      Nr   = r_d.getSize () - 1;
      rr   = new double[Nr+1];
      PsiN = new double[Nr+1];
      Psi  = new double[Nr+1];
      f    = new double[Nr+1];
      g2   = new double[Nr+1];
      p2   = new double[Nr+1];
      pp   = new double[Nr+1];
      ppp  = new double[Nr+1];
      q    = new double[Nr+1];
      s    = new double[Nr+1];
      s2   = new double[Nr+1];
      s0   = new double[Nr+1];
      S1   = new double[Nr+1];
      S3   = new double[Nr+1];
      S4   = new double[Nr+1];
      P1   = new double[Nr+1];
      P2   = new double[Nr+1];
      P1a  = new double[Nr+1];
      P2a  = new double[Nr+1];
      P3   = new double[Nr+1];
      P4   = new double[Nr+1];

      r_x.  getVar (rr);
      P_x.  getVar (PsiN);
      PP_x. getVar (Psi);
      f_x.  getVar (f);
      g2_x. getVar (g2);
      p2_x. getVar (p2);
      pp_x. getVar (pp);
      ppp_x.getVar (ppp);
      q_x  .getVar (q);
      s_x  .getVar (s);
      s2_x .getVar (s2);
      s0_x .getVar (s0);
      S1_x .getVar (S1);
      S3_x .getVar (S3);
      S4_x .getVar (S4);
      P1_x .getVar (P1);
      P2_x .getVar (P2);
      P1a_x.getVar (P1a);
      P2a_x.getVar (P2a);
      P3_x .getVar (P3);
      P4_x .getVar (P4);

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
        
      delete[] para;
      delete[] Hndata; delete[] Hnpdata; delete[] Vndata; delete[] Vnpdata;
      delete[] RRdata; delete[] ZZdata;  delete[] rrdata; delete[] ttdata;
    }
  catch (NcException& e)
    {
      printf ("Error reading data from netcdf file Outputs/Equilibrium/Equilbrium.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }
}

// ###############################################
// Function to write stability data to netcdf file
// ###############################################
void TJ::WriteNetcdf ()
{
  printf ("Writing stability data to netcdf file Outputs/TJ/TJ.nc:\n");
 
  double Input[23];

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
  double* km      = new double[(Nr+1)*J];

  double* Pvac_r  = new double[J*J];
  double* Pvac_i  = new double[J*J];
  double* Rvac_r  = new double[J*J];
  double* Rvac_i  = new double[J*J];
  double* Amat_r  = new double[J*J];
  double* Amat_i  = new double[J*J];
  double* Aant_r  = new double[J*J];
  double* Aant_i  = new double[J*J];
  double* Hmat_r  = new double[J*J];
  double* Hmat_i  = new double[J*J];

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
  double* PPV_r   = new double[nres*Nf*(Nw+1)];
  double* PPV_i   = new double[nres*Nf*(Nw+1)];
  double* ZZV_r   = new double[nres*Nf*(Nw+1)];
  double* ZZV_i   = new double[nres*Nf*(Nw+1)];
  double* Fmat_r  = new double[nres*nres];
  double* Fmat_i  = new double[nres*nres];
  double* Emat_r  = new double[nres*nres];
  double* Emat_i  = new double[nres*nres];
  double* Eant_r  = new double[nres*nres];
  double* Eant_i  = new double[nres*nres];
  double* Fvec_r  = new double[nres*nres];
  double* Fvec_i  = new double[nres*nres];
  double* Psix_r  = new double[J];
  double* Psix_i  = new double[J];

  double* Xi_r    = new double[J];
  double* Xi_i    = new double[J];
  double* Up_r    = new double[J];
  double* Up_i    = new double[J];
  double* Chi_r   = new double[nres];
  double* Chi_i   = new double[nres];
  double* Chi_a   = new double[nres];
  double* Delta_r = new double[nres];
  double* PPR_r   = new double[J*NDIAG];
  double* PPR_i   = new double[J*NDIAG];
  double* ZZR_r   = new double[J*NDIAG];
  double* ZZR_i   = new double[J*NDIAG];
  double* PPRV_r  = new double[Nf*(Nw+1)];
  double* PPRV_i  = new double[Nf*(Nw+1)];
  double* ZZRV_r  = new double[Nf*(Nw+1)];
  double* ZZRV_i  = new double[Nf*(Nw+1)];

  double* Psii_r  = new double[J*J*NDIAG];
  double* Psii_i  = new double[J*J*NDIAG];
  double* Zi_r    = new double[J*J*NDIAG];
  double* Zi_i    = new double[J*J*NDIAG];
  double* Xii_r   = new double[J*J*NDIAG];
  double* Xii_i   = new double[J*J*NDIAG];
  double* Chii_r  = new double[J*J*NDIAG];
  double* Chii_i  = new double[J*J*NDIAG];
  double* Umat_r  = new double[J*J];
  double* Umat_i  = new double[J*J];
  double* Uant_r  = new double[J*J];
  double* Uant_i  = new double[J*J];
  double* Uvec_r  = new double[J*J];
  double* Uvec_i  = new double[J*J];
  double* U1mat_r = new double[J*J];
  double* U1mat_i = new double[J*J];
  double* Uvec1_r = new double[J*J];
  double* Uvec1_i = new double[J*J];
  double* Psie_r  = new double[J*J*NDIAG];
  double* Psie_i  = new double[J*J*NDIAG];
  double* Ze_r    = new double[J*J*NDIAG];
  double* Ze_i    = new double[J*J*NDIAG];
  double* Xie_r   = new double[J*J*NDIAG];
  double* Xie_i   = new double[J*J*NDIAG];
  double* xvals   = new double[J*NDIAG];

  double* Psiy_r  = new double[J*(Nw+1)];
  double* Psiy_i  = new double[J*(Nw+1)];
  double* Jy_r    = new double[J*(Nw+1)];
  double* Jy_i    = new double[J*(Nw+1)];
  double* Xiy_r   = new double[J*(Nw+1)];
  double* Xiy_i   = new double[J*(Nw+1)];

  double* Psis_r  = new double[Nw+1];
  double* Psis_i  = new double[Nw+1];
  double* Psiz_r  = new double[Nw+1];
  double* Psiz_i  = new double[Nw+1];

  double* gammax_r = new double[J];
  double* gammax_i = new double[J];
  double* gamma_r  = new double[J];
  double* gamma_i  = new double[J];

  Input[0]  = double (NTOR);
  Input[1]  = double (MMIN);
  Input[2]  = double (MMAX);
  Input[3]  = EPS;
  Input[4]  = DEL;
  Input[5]  = double (NFIX);
  Input[6]  = double (NDIAG);
  Input[7]  = NULC;
  Input[8]  = double (ITERMAX);
  Input[9]  = double (FREE);
  Input[10] = acc;
  Input[11] = h0;
  Input[12] = hmin;
  Input[13] = hmax;
  Input[14] = EPSF;
  Input[15] = B0;
  Input[16] = R0;
  Input[17] = n0;
  Input[18] = alpha;
  Input[19] = Zeff;
  Input[20] = Mion;
  Input[21] = Chip;
  Input[22] = Teped;
  
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
  for (int i = 0; i <= Nr; i++)
    for (int j = 0; j < J; j++)
      {
	km[cnt] = Getkm (rr[i], MPOL[j]);
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
	Amat_r[cnt] = Amat(j, jp).real();
	Amat_i[cnt] = Amat(j, jp).imag();
	Aant_r[cnt] = Aant(j, jp).real();
	Aant_i[cnt] = Aant(j, jp).imag();
	Hmat_r[cnt] = Hmat(j, jp).real();
	Hmat_i[cnt] = Hmat(j, jp).imag();
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

  cnt = 0;
  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Fmat_r[cnt] = real (Fmat(j, jp));
	Fmat_i[cnt] = imag (Fmat(j, jp));
	Emat_r[cnt] = real (Emat(j, jp));
	Emat_i[cnt] = imag (Emat(j, jp));
	Eant_r[cnt] = 0.5 * real (Emat(j, jp) - conj (Emat(jp ,j)));
	Eant_i[cnt] = 0.5 * imag (Emat(j, jp) - conj (Emat(jp ,j)));
	Fvec_r[cnt] = real (Fvec(j, jp));
	Fvec_i[cnt] = imag (Fvec(j, jp));
	cnt++;
      }

  for (int j = 0; j < J; j++)
    {
      Psix_r[j] = real (Psix[j]);
      Psix_i[j] = imag (Psix[j]);
      Xi_r  [j] = real (Xi[j]);
      Xi_i  [j] = imag (Xi[j]);
      Up_r  [j] = real (Upsilon[j]);
      Up_i  [j] = imag (Upsilon[j]);
    }

  for (int j = 0; j < nres; j++)
    {
      Chi_r[j] = real (Chi[j]);
      Chi_i[j] = imag (Chi[j]);
      Chi_a[j] = abs (Chi[j]);
    }

  for (int j = 0; j < nres; j++)
    Delta_r[j] = real (Emat(j, j));

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int i = 0; i < NDIAG; i++)
      {
	PPR_r[cnt] = real (Psirmp(j, i));
	PPR_i[cnt] = imag (Psirmp(j, i));
	ZZR_r[cnt] = real (Zrmp  (j, i));
	ZZR_i[cnt] = imag (Zrmp  (j, i));
	cnt++;
      }

  cnt = 0;
  for (int i = 0; i < Nf; i++)
    for (int l = 0; l <= Nw; l++)
      {
	PPRV_r[cnt] = real (Psirv(i, l));
	PPRV_i[cnt] = imag (Psirv(i, l));
	ZZRV_r[cnt] = real (Zrv  (i, l));
	ZZRV_i[cnt] = imag (Zrv  (i, l));
	cnt++;
      }

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
      {
	Psii_r [cnt] = real (Psii (j, jp, i));
	Psii_i [cnt] = imag (Psii (j, jp, i));
	Zi_r   [cnt] = real (Zi   (j, jp, i));
	Zi_i   [cnt] = imag (Zi   (j, jp, i));
	Xii_r  [cnt] = real (Xii  (j, jp, i));
	Xii_i  [cnt] = imag (Xii  (j, jp, i));
	Chii_r [cnt] = real (Chii (j, jp, i));
	Chii_i [cnt] = imag (Chii (j, jp, i));
	cnt++;
      }
  
  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Umat_r[cnt]  = Umat(j, jp) .real();
	Umat_i[cnt]  = Umat(j, jp) .imag();
	Uant_r[cnt]  = Uant(j, jp) .real();
	Uant_i[cnt]  = Uant(j, jp) .imag();
	Uvec_r[cnt]  = Uvec(j, jp) .real();
	Uvec_i[cnt]  = Uvec(j, jp) .imag();
	U1mat_r[cnt] = U1mat(j, jp).real();
	U1mat_i[cnt] = U1mat(j, jp).imag();
	Uvec1_r[cnt] = U1vec(j, jp).real();
	Uvec1_i[cnt] = U1vec(j, jp).imag(); 
	cnt++;
      }
  
  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  Psie_r[cnt] = real (Psie(j, jp, i));
	  Psie_i[cnt] = imag (Psie(j, jp, i));
	  Ze_r  [cnt] = real (Ze  (j, jp, i));
	  Ze_i  [cnt] = imag (Ze  (j, jp, i));
	  Xie_r [cnt] = real (Xie (j, jp, i));
	  Xie_i [cnt] = imag (Xie (j, jp, i));
	  cnt++;
	}

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int i = 0; i < NDIAG; i++)
      {
	xvals[cnt] = lvals(j, i);
	cnt++;
      }
  
  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int i = 0; i <= Nw; i++)
      {
	Psiy_r[cnt] = real (Psiy(j, i));
	Psiy_i[cnt] = imag (Psiy(j, i));
	Jy_r  [cnt] = real (Jy  (j, i));
	Jy_i  [cnt] = imag (Jy  (j, i));
	Xiy_r [cnt] = real (Xiy (j, i));
	Xiy_i [cnt] = imag (Xiy (j, i));
	cnt++;
      }

  for (int i = 0; i <= Nw; i++)
    {
      Psis_r[i] = real (Psirmps[i]);
      Psis_i[i] = imag (Psirmps[i]);
      Psiz_r[i] = real (Psixs  [i]);
      Psiz_i[i] = imag (Psixs  [i]);
    }

  for (int j = 0; j < J; j++)
    {
      gammax_r[j] = real (gammax[j]);
      gammax_i[j] = imag (gammax[j]);
      gamma_r [j] = real (gamma [j]);
      gamma_i [j] = imag (gamma [j]);
    }
    
  try
    {
      NcFile dataFile ("../Outputs/TJ/TJ.nc", NcFile::replace);

      NcDim i_d = dataFile.addDim ("Ni",    23);
      NcDim r_d = dataFile.addDim ("Nr",    Nr+1);
      NcDim s_d = dataFile.addDim ("Ns",    Ns+1);
      NcDim x_d = dataFile.addDim ("nres",  nres);
      NcDim j_d = dataFile.addDim ("J",     J);
      NcDim k_d = dataFile.addDim ("K",     K);
      NcDim d_d = dataFile.addDim ("ndiag", NDIAG);
      NcDim f_d = dataFile.addDim ("Nf",    Nf);
      NcDim w_d = dataFile.addDim ("Nw",    Nw+1);
      NcDim c_d = dataFile.addDim ("Ncoil", ncoil);

      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> matrix_d;
      matrix_d.push_back (r_d);
      matrix_d.push_back (j_d);
      matrix_d.push_back (j_d);

      vector<NcDim> km_d;
      km_d.push_back (r_d);
      km_d.push_back (j_d);

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

      vector<NcDim> e_d;
      e_d.push_back (x_d);
      e_d.push_back (x_d);

      vector<NcDim> rmp_d;
      rmp_d.push_back (j_d);
      rmp_d.push_back (d_d);

      vector<NcDim> vr_d;
      vr_d.push_back (f_d);
      vr_d.push_back (w_d);

      vector<NcDim> ideal_d;
      ideal_d.push_back (j_d);
      ideal_d.push_back (j_d);
      ideal_d.push_back (d_d);

      vector<NcDim> surface_d;
      surface_d.push_back (j_d);
      surface_d.push_back (w_d);

      NcVar i_x = dataFile.addVar ("InputParameters", ncDouble, i_d);
      i_x.putVar (Input);
      
      NcVar r_x   = dataFile.addVar ("r",   ncDouble, r_d);
      r_x.putVar (rr);
      NcVar pp_x  = dataFile.addVar ("pp",  ncDouble, r_d);
      pp_x.putVar (pp);
      NcVar ppp_x = dataFile.addVar ("ppp", ncDouble, r_d);
      ppp_x.putVar (ppp);
      NcVar q_x   = dataFile.addVar ("q",   ncDouble, r_d);
      q_x.putVar (q);
      NcVar s_x   = dataFile.addVar ("s",   ncDouble, r_d);
      s_x.putVar (s);
      NcVar s2_x  = dataFile.addVar ("s2",  ncDouble, r_d);
      s2_x.putVar (s2);
      NcVar S1_x  = dataFile.addVar ("S1",  ncDouble, r_d);
      S1_x.putVar (S1);
      NcVar P1_x  = dataFile.addVar ("P1",  ncDouble, r_d);
      P1_x.putVar (P1);
      NcVar P2_x  = dataFile.addVar ("P2",  ncDouble, r_d);
      P2_x.putVar (P2);
      NcVar P3_x  = dataFile.addVar ("P3",  ncDouble, r_d);
      P3_x.putVar (P3);
 
      NcVar Hn_x  = dataFile.addVar ("Hn",  ncDouble, shape_d);
      Hn_x.putVar (HHfunc.data());
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (HPfunc.data());
      NcVar Vn_x  = dataFile.addVar ("Vn",  ncDouble, shape_d);
      Vn_x.putVar (VVfunc.data());
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (VPfunc.data());

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

      NcVar ltest_x  = dataFile.addVar ("Ltest",  ncDouble, matrix_d);
      ltest_x.putVar (Ltest);
      NcVar mntest_x = dataFile.addVar ("MNtest", ncDouble, matrix_d);
      mntest_x.putVar (MNtest);
      NcVar ptest_x  = dataFile.addVar ("Ptest",  ncDouble, matrix_d);
      ptest_x.putVar (Ptest);
      NcVar km_x     = dataFile.addVar ("km",     ncDouble, km_d);
      km_x.putVar (km);

      NcVar pvacr_x = dataFile.addVar ("Pvac_r", ncDouble, vacuum_d);
      pvacr_x.putVar (Pvac_r);
      NcVar pvaci_x = dataFile.addVar ("Pvac_i", ncDouble, vacuum_d);
      pvaci_x.putVar (Pvac_i);
      NcVar rvacr_x = dataFile.addVar ("Rvac_r", ncDouble, vacuum_d);
      rvacr_x.putVar (Rvac_r);
      NcVar rvaci_x = dataFile.addVar ("Rvac_i", ncDouble, vacuum_d);
      rvaci_x.putVar (Rvac_i);
      
      NcVar avacr_x = dataFile.addVar ("Amat_r", ncDouble, vacuum_d);
      avacr_x.putVar (Amat_r);
      NcVar avaci_x = dataFile.addVar ("Amat_i", ncDouble, vacuum_d);
      avaci_x.putVar (Amat_i);
      NcVar aantr_x = dataFile.addVar ("Aant_r", ncDouble, vacuum_d);
      aantr_x.putVar (Aant_r);
      NcVar aanti_x = dataFile.addVar ("Aant_i", ncDouble, vacuum_d);
      aanti_x.putVar (Aant_i);
 
      NcVar hmatr_x = dataFile.addVar ("Hmat_r", ncDouble, vacuum_d);
      hmatr_x.putVar (Hmat_r);
      NcVar hmati_x = dataFile.addVar ("Hmat_i", ncDouble, vacuum_d);
      hmati_x.putVar (Hmat_i);
   
      NcVar rgrid_x = dataFile.addVar ("r_grid",    ncDouble, d_d);
      rgrid_x.putVar (Rgrid);
      NcVar pgrid_x = dataFile.addVar ("PsiN_grid", ncDouble, d_d);
      pgrid_x.putVar (Pgrid);
      NcVar hode_x  = dataFile.addVar ("h_ode",     ncDouble, d_d);
      hode_x.putVar (hode);
      NcVar eode_x  = dataFile.addVar ("err_ode",   ncDouble, d_d);
      eode_x.putVar (eode);

      NcVar t_x      = dataFile.addVar  ("theta",  ncDouble, w_d);
      t_x.putVar (tbound);
      NcVar cmu_x     = dataFile.addVar ("cosmu",  ncDouble, w_d);
      cmu_x.putVar (cmu);
      NcVar e_x      = dataFile.addVar  ("eta",    ncDouble, w_d);
      e_x.putVar (eeta);
      NcVar ceta_x   = dataFile.addVar  ("coseta", ncDouble, w_d);
      ceta_x.putVar (ceta);
      NcVar seta_x   = dataFile.addVar  ("sineta", ncDouble, w_d);
      seta_x.putVar (seta);
      NcVar R2grgz_x = dataFile.addVar  ("R2grgz", ncDouble, w_d);
      R2grgz_x.putVar (R2grgz);
      NcVar R2grge_x = dataFile.addVar  ("R2grge", ncDouble, w_d);
      R2grge_x.putVar (R2grge);

      NcVar ttest_x = dataFile.addVar ("Torque_test", ncDouble, torque_d);
      ttest_x.putVar (Ttest.data());
      NcVar pnorm_x = dataFile.addVar ("Psi_norm",    ncDouble, torque_d);
      pnorm_x.putVar (Pnorm.data());
      NcVar znorm_x = dataFile.addVar ("Z_norm",      ncDouble, torque_d);
      znorm_x.putVar (Znorm.data());

      NcVar pppsir_x = dataFile.addVar ("Psi_r", ncDouble, soln_d);
      pppsir_x.putVar (PPPsi_r);
      NcVar pppsii_x = dataFile.addVar ("Psi_i", ncDouble, soln_d);
      pppsii_x.putVar (PPPsi_i);
      NcVar zzzr_x   = dataFile.addVar ("Z_r",   ncDouble, soln_d);
      zzzr_x.putVar (ZZZ_r);
      NcVar zzzi_x   = dataFile.addVar ("Z_i",   ncDouble, soln_d);
      zzzi_x.putVar (ZZZ_i);

      NcVar ppfr_x = dataFile.addVar ("Psi_full_r", ncDouble, full_d);
      ppfr_x.putVar (PPF_r);
      NcVar ppfi_x = dataFile.addVar ("Psi_full_i", ncDouble, full_d);
      ppfi_x.putVar (PPF_i);
      NcVar zzfr_x = dataFile.addVar ("Z_full_r",   ncDouble, full_d);
      zzfr_x.putVar (ZZF_r);
      NcVar zzfi_x = dataFile.addVar ("Z_full_i",   ncDouble, full_d);
      zzfi_x.putVar (ZZF_i);

      NcVar ppur_x = dataFile.addVar ("Psi_unrc_r", ncDouble, full_d);
      ppur_x.putVar (PPU_r);
      NcVar ppui_x = dataFile.addVar ("Psi_unrc_i", ncDouble, full_d);
      ppui_x.putVar (PPU_i);
      NcVar zzur_x = dataFile.addVar ("Z_unrc_r",   ncDouble, full_d);
      zzur_x.putVar (ZZU_r);
      NcVar zzui_x = dataFile.addVar ("Z_unrc_i",   ncDouble, full_d);
      zzui_x.putVar (ZZU_i);

      NcVar mres_x = dataFile.addVar   ("m_res",  ncInt,    x_d);
      mres_x.putVar (mres);
      NcVar rres_x = dataFile.addVar   ("r_res",  ncDouble, x_d);
      rres_x.putVar (rres);
      NcVar sres_x = dataFile.addVar   ("s_res",  ncDouble, x_d);
      sres_x.putVar (sres);
      NcVar dires_x = dataFile.addVar  ("DI_res", ncDouble, x_d);
      dires_x.putVar (DIres);
      NcVar drres_x = dataFile.addVar  ("DR_res", ncDouble, x_d);
      drres_x.putVar (DRres);
      NcVar flarge_x = dataFile.addVar ("Flarge", ncDouble, x_d);
      flarge_x.putVar (Flarge);
      NcVar fsmall_x = dataFile.addVar ("Fsmall", ncDouble, x_d);
      fsmall_x.putVar (Fsmall);

      NcVar tf_x = dataFile.addVar ("Torque_full", ncDouble, t_d);
      tf_x.putVar (Tf.data());
      NcVar tu_x = dataFile.addVar ("Torque_unrc", ncDouble, t_d);
      tu_x.putVar (Tu.data());

      NcVar tfull_x = dataFile.addVar ("Torque_pair_full", ncDouble, tt_d);
      tfull_x.putVar (Tfull.data());
      NcVar tunrc_x = dataFile.addVar ("Torque_pair_unrc", ncDouble, tt_d);
      tunrc_x.putVar (Tunrc.data());

      NcVar ppvr_x = dataFile.addVar ("Psi_unrc_eig_r", ncDouble, v_d);
      ppvr_x.putVar (PPV_r);
      NcVar ppvi_x = dataFile.addVar ("Psi_unrc_eig_i", ncDouble, v_d);
      ppvi_x.putVar (PPV_i);
      NcVar zzvr_x = dataFile.addVar ("Z_unrc_eig_r",   ncDouble, v_d);
      zzvr_x.putVar (ZZV_r);
      NcVar zzvi_x = dataFile.addVar ("Z_unrc_eig_i",   ncDouble, v_d);
      zzvi_x.putVar (ZZV_i);

      NcVar fmatr_x = dataFile.addVar ("Fmat_r", ncDouble, e_d);
      fmatr_x.putVar (Fmat_r);
      NcVar fmati_x = dataFile.addVar ("Fmat_i", ncDouble, e_d);
      fmati_x.putVar (Fmat_i);
      NcVar ematr_x = dataFile.addVar ("Emat_r", ncDouble, e_d);
      ematr_x.putVar (Emat_r);
      NcVar emati_x = dataFile.addVar ("Emat_i", ncDouble, e_d);
      emati_x.putVar (Emat_i);
      NcVar eantr_x = dataFile.addVar ("Eant_r", ncDouble, e_d);
      eantr_x.putVar (Eant_r);
      NcVar eanti_x = dataFile.addVar ("Eant_i", ncDouble, e_d);
      eanti_x.putVar (Eant_i);
      NcVar Fval_x  = dataFile.addVar ("Fval",   ncDouble, x_d);
      Fval_x.putVar (Fval);
      NcVar Fvecr_x = dataFile.addVar ("Fvec_r", ncDouble, e_d);
      Fvecr_x.putVar (Fvec_r);
      NcVar Fveci_x = dataFile.addVar ("Fvec_i", ncDouble, e_d);
      Fveci_x.putVar (Fvec_i);

      NcVar rcoil_x = dataFile.addVar ("Rcoil",     ncDouble, c_d);
      rcoil_x.putVar (Rcoil);
      NcVar zcoil_x = dataFile.addVar ("Zcoil",     ncDouble, c_d);
      zcoil_x.putVar (Zcoil);
      NcVar icoil_x = dataFile.addVar ("Icoil",     ncDouble, c_d);
      icoil_x.putVar (Icoil);
      NcVar psix_r  = dataFile.addVar ("Psi^x_r",   ncDouble, j_d);
      psix_r.putVar (Psix_r);
      NcVar psix_i  = dataFile.addVar ("Psi^x_i",   ncDouble, j_d);
      psix_i.putVar (Psix_i);
      NcVar xi_r    = dataFile.addVar ("Xi_r",      ncDouble, j_d);
      xi_r.putVar (Xi_r);
      NcVar xi_i    = dataFile.addVar ("Xi_i",      ncDouble, j_d);
      xi_i.putVar (Xi_i);
      NcVar up_r    = dataFile.addVar ("Upsilon_r", ncDouble, j_d);
      up_r.putVar (Up_r);
      NcVar up_i    = dataFile.addVar ("Upsilon_i", ncDouble, j_d);
      up_i.putVar (Up_i);
      NcVar chi_r   = dataFile.addVar ("Chi_r",     ncDouble, x_d);
      chi_r.putVar (Chi_r);
      NcVar chi_i   = dataFile.addVar ("Chi_i",     ncDouble, x_d);
      chi_i.putVar (Chi_i);
      NcVar chi_a   = dataFile.addVar ("Chi_a",     ncDouble, x_d);
      chi_a.putVar (Chi_a);

      NcVar Te_r  = dataFile.addVar ("Te",         ncDouble, x_d);
      Te_r.putVar (Teres);
      NcVar S13_r = dataFile.addVar ("S13",        ncDouble, x_d);
      S13_r.putVar (S13res);
      NcVar tau_r = dataFile.addVar ("tau",        ncDouble, x_d);
      tau_r.putVar (taures);
      NcVar ie_r  = dataFile.addVar ("iota_e",     ncDouble, x_d);
      ie_r.putVar (ieres);
      NcVar Qe_r  = dataFile.addVar ("Qe",         ncDouble, x_d);
      Qe_r.putVar (Qeres);
      NcVar Qi_r  = dataFile.addVar ("Qi",         ncDouble, x_d);
      Qi_r.putVar (Qires);
      NcVar QE_r  = dataFile.addVar ("QE",         ncDouble, x_d);
      QE_r.putVar (QEres);
      NcVar D_r   = dataFile.addVar ("D",          ncDouble, x_d);
      D_r.putVar (Dres);
      NcVar Pm_r  = dataFile.addVar ("Pphi",       ncDouble, x_d);
      Pm_r.putVar (Pmres);
      NcVar Pe_r  = dataFile.addVar ("Pperp",      ncDouble, x_d);
      Pe_r.putVar (Peres);
      NcVar Dc_r  = dataFile.addVar ("Delta_crit", ncDouble, x_d);
      Dc_r.putVar (Dcres);
      NcVar De_r  = dataFile.addVar ("Delta",      ncDouble, x_d);
      De_r.putVar (Delta_r);

      NcVar pprr_x = dataFile.addVar ("Psi_rmp_r", ncDouble, rmp_d);
      pprr_x.putVar (PPR_r);
      NcVar ppri_x = dataFile.addVar ("Psi_rmp_i", ncDouble, rmp_d);
      ppri_x.putVar (PPR_i);
      NcVar zzrr_x = dataFile.addVar ("Z_rmp_r",   ncDouble, rmp_d);
      zzrr_x.putVar (ZZR_r);
      NcVar zzri_x = dataFile.addVar ("Z_rmp_i",   ncDouble, rmp_d);
      zzri_x.putVar (ZZR_i);

      NcVar pprvr_x = dataFile.addVar ("Psi_rmp_eig_r", ncDouble, vr_d);
      pprvr_x.putVar (PPRV_r);
      NcVar pprvi_x = dataFile.addVar ("Psi_rmp_eig_i", ncDouble, vr_d);
      pprvi_x.putVar (PPRV_i);
      NcVar zzrvr_x = dataFile.addVar ("Z_rmp_eig_r",   ncDouble, vr_d);
      zzrvr_x.putVar (ZZRV_r);
      NcVar zzrvi_x = dataFile.addVar ("Z_rmp_eig_i",   ncDouble, vr_d);
      zzrvi_x.putVar (ZZRV_i);

      NcVar psiir_x = dataFile.addVar ("Psi_i_r", ncDouble, ideal_d);
      psiir_x.putVar (Psii_r);
      NcVar psiii_x = dataFile.addVar ("Psi_i_i", ncDouble, ideal_d);
      psiii_x.putVar (Psii_i);
      NcVar zir_x   = dataFile.addVar ("Z_i_r",   ncDouble, ideal_d);
      zir_x.putVar (Zi_r);
      NcVar zii_x   = dataFile.addVar ("Z_i_i",   ncDouble, ideal_d);
      zii_x.putVar (Zi_i);
      NcVar xir_x   = dataFile.addVar ("Xi_i_r",  ncDouble, ideal_d);
      xir_x.putVar (Xii_r);
      NcVar xii_x   = dataFile.addVar ("Xi_i_i",  ncDouble, ideal_d);
      xii_x.putVar (Xii_i);
      NcVar chiir_x = dataFile.addVar ("Chi_i_r", ncDouble, ideal_d);
      chiir_x.putVar (Chii_r);
      NcVar chiii_x = dataFile.addVar ("Chi_i_i", ncDouble, ideal_d);
      chiii_x.putVar (Chii_i);
      NcVar lvals_x = dataFile.addVar ("lambda",  ncDouble, rmp_d);
      lvals_x.putVar (xvals);
      
      NcVar wvacr_x  = dataFile.addVar ("Umat_r",  ncDouble, vacuum_d);
      wvacr_x.putVar (Umat_r);
      NcVar wvaci_x  = dataFile.addVar ("Umat_i",  ncDouble, vacuum_d);
      wvaci_x.putVar (Umat_i);
      NcVar wantr_x  = dataFile.addVar ("Uant_r",  ncDouble, vacuum_d);
      wantr_x.putVar (Uant_r);
      NcVar wanti_x  = dataFile.addVar ("Uant_i",  ncDouble, vacuum_d);
      wanti_x.putVar (Uant_i);
      NcVar wval_x   = dataFile.addVar ("Wval",    ncDouble, j_d);
      wval_x.putVar  (Wval);
      NcVar vval_x   = dataFile.addVar ("Vval",    ncDouble, j_d);
      vval_x.putVar  (Vval);
      NcVar uval_x   = dataFile.addVar ("Uval",    ncDouble, j_d);
      uval_x.putVar  (Uval);
      NcVar wvecr_x  = dataFile.addVar ("Uvec_r",  ncDouble, vacuum_d);
      wvecr_x.putVar (Uvec_r);
      NcVar w1vacr_x = dataFile.addVar ("U1mat_r", ncDouble, vacuum_d);
      w1vacr_x.putVar (U1mat_r);
      NcVar w1vaci_x = dataFile.addVar ("U1mat_i", ncDouble, vacuum_d);
      w1vaci_x.putVar (U1mat_i);
      NcVar wveci_x  = dataFile.addVar ("Uvec_i",  ncDouble, vacuum_d);
      wveci_x.putVar (Uvec_i);
      NcVar wvec1r_x = dataFile.addVar ("Uvec1_r", ncDouble, vacuum_d);
      wvec1r_x.putVar (Uvec1_r);
      NcVar wvec1i_x = dataFile.addVar ("Uvec1_i", ncDouble, vacuum_d);
      wvec1i_x.putVar (Uvec1_i);

      NcVar psier_x = dataFile.addVar ("Psi_e_r", ncDouble, ideal_d);
      psier_x.putVar (Psie_r);
      NcVar psiei_x = dataFile.addVar ("Psi_e_i", ncDouble, ideal_d);
      psiei_x.putVar (Psie_i);
      NcVar zer_x   = dataFile.addVar ("Z_e_r",   ncDouble, ideal_d);
      zer_x.putVar (Ze_r);
      NcVar zei_x   = dataFile.addVar ("Z_e_i",   ncDouble, ideal_d);
      zei_x.putVar (Ze_i);
      NcVar xer_x   = dataFile.addVar ("Xi_e_r",  ncDouble, ideal_d);
      xer_x.putVar (Xie_r);
      NcVar xei_x   = dataFile.addVar ("Xi_e_i",  ncDouble, ideal_d);
      xei_x.putVar (Xie_i);

      NcVar dW_x  = dataFile.addVar ("delta_W",   ncDouble, j_d);
      dW_x.putVar (deltaW);
      NcVar dWv_x = dataFile.addVar ("delta_W_p", ncDouble, j_d);
      dWv_x.putVar (deltaWp);
      NcVar dWp_x = dataFile.addVar ("delta_W_v", ncDouble, j_d);
      dWp_x.putVar (deltaWv);

      NcVar psiyr_x = dataFile.addVar ("Psi_surface_r", ncDouble, surface_d);
      psiyr_x.putVar (Psiy_r);
      NcVar psiyi_x = dataFile.addVar ("Psi_surface_i", ncDouble, surface_d);
      psiyi_x.putVar (Psiy_i);
      NcVar jyr_x   = dataFile.addVar ("J_surface_r",   ncDouble, surface_d);
      jyr_x.putVar (Jy_r);
      NcVar jyi_x   = dataFile.addVar ("J_surface_i",   ncDouble, surface_d);
      jyi_x.putVar (Jy_i);
      NcVar xyr_x   = dataFile.addVar ("Xi_surface_r",  ncDouble, surface_d);
      xyr_x.putVar (Xiy_r);
      NcVar xyi_x   = dataFile.addVar ("Xi_surface_i",  ncDouble, surface_d);
      xyi_x.putVar (Xiy_i);
 
      NcVar psisr_x = dataFile.addVar ("Psi_rmp_surface_r", ncDouble, w_d);
      psisr_x.putVar (Psis_r);
      NcVar psisi_x = dataFile.addVar ("Psi_rmp_surface_i", ncDouble, w_d);
      psisi_x.putVar (Psis_i);
      NcVar psixr_x = dataFile.addVar ("Psi_x_surface_r",   ncDouble, w_d);
      psixr_x.putVar (Psiz_r);
      NcVar psixi_x = dataFile.addVar ("Psi_x_surface_i",   ncDouble, w_d);
      psixi_x.putVar (Psiz_i);

      NcVar gammaxr_x = dataFile.addVar ("gammax_r", ncDouble, j_d);
      gammaxr_x.putVar (gammax_r);
      NcVar gammaxi_x = dataFile.addVar ("gammax_i", ncDouble, j_d);
      gammaxi_x.putVar (gammax_i);
      NcVar gammar_x = dataFile.addVar  ("gamma_r",  ncDouble, j_d);
      gammar_x.putVar (gamma_r);
      NcVar gammai_x = dataFile.addVar  ("gamma_i",  ncDouble, j_d);
      gammai_x.putVar (gamma_i);
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/TJ/TJ.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] Lmmp_r; delete[] Mmmp_r; delete[] Nmmp_r; delete[] Pmmp_r;
  delete[] Lmmp_i; delete[] Mmmp_i; delete[] Nmmp_i; delete[] Pmmp_i;
  delete[] Ltest;  delete[] MNtest; delete[] Ptest;  delete[] km;

  delete[] Pvac_r; delete[] Pvac_i; delete[] Rvac_r; delete[] Rvac_i;
  delete[] Amat_r; delete[] Amat_i; delete[] Aant_r; delete[] Aant_i;
  delete[] Hmat_r; delete[] Hmat_i; 

  delete[] PPPsi_r; delete[] PPPsi_i; delete[] ZZZ_r;   delete[] ZZZ_i;
  delete[] PPF_r;   delete[] PPF_i;   delete[] ZZF_r;   delete[] ZZF_i;
  delete[] PPU_r;   delete[] PPU_i;   delete[] ZZU_r;   delete[] ZZU_i;

  delete[] PPV_r;   delete[] PPV_i;  delete[] ZZV_r;  delete[] ZZV_i;
  delete[] Emat_r;  delete[] Emat_i; delete[] Eant_r; delete[] Eant_i;
  delete[] Psix_r;  delete[] Psix_i; delete[] Xi_r;   delete[] Xi_i; 
  delete[] Up_r;    delete[] Up_i;   delete[] Chi_r;  delete[] Chi_i;
  delete[] Delta_r; delete[] Fvec_r; delete[] Fvec_i; delete[] Fmat_r;
  delete[] Fmat_i;

  delete[] PPR_r;  delete[] PPR_i;  delete[] ZZR_r;  delete[] ZZR_i;
  delete[] PPRV_r; delete[] PPRV_i; delete[] ZZRV_r; delete[] ZZRV_i;

  delete[] Psii_r;   delete[] Psii_i;   delete[] Zi_r;    delete[] Zi_i;
  delete[] Xii_r;    delete[] Xii_i;    delete[] Uvec_r;  delete[] Uvec_i;
  delete[] Umat_r;   delete[] Umat_i;   delete[] Uant_r;  delete[] Uant_i;
  delete[] Psie_r;   delete[] Psie_i;   delete[] Ze_r;    delete[] Ze_i;
  delete[] Xie_r;    delete[] Xie_i;    delete[] Psiy_r;  delete[] Psiy_i;
  delete[] Jy_r;     delete[] Jy_i;     delete[] Xiy_r;   delete[] Xiy_i;
  delete[] Psis_r;   delete[] Psis_i;   delete[] Psiz_r;  delete[] Psiz_i;
  delete[] gammax_r; delete[] gammax_i; delete[] gamma_r; delete[] gamma_i;
  delete[] Uvec1_r;  delete[] Uvec1_i;  delete[] U1mat_r; delete[] U1mat_i;
  delete[] Chii_r;   delete[] Chii_i;   
  delete[] xvals;
}
