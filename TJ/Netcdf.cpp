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
      bw   = para[2];
      
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
      NcVar S2_x  = dataFile.getVar ("S2");
      NcVar S3_x  = dataFile.getVar ("S3");
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P3_x  = dataFile.getVar ("P3");
      NcDim r_d   = r_x.getDim (0);

      Nr   = r_d.getSize () - 1;
      Nint = Nr;
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
      S2   = new double[Nr+1];
      S3   = new double[Nr+1];
      P1   = new double[Nr+1];
      P2   = new double[Nr+1];
      P3   = new double[Nr+1];

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
      S2_x .getVar (S2);
      S3_x .getVar (S3);
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

      NcVar tbound_x = dataFile.getVar ("tbound");
      NcVar Rbound_x = dataFile.getVar ("Rbound");
      NcVar Zbound_x = dataFile.getVar ("Zbound");
      NcVar dRdthe_x = dataFile.getVar ("dRdtheta");
      NcVar dZdthe_x = dataFile.getVar ("dZdtheta");
      NcDim w_d      = tbound_x.getDim (0);

      Nw = w_d.getSize () - 1;
 
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

      NcVar wwall_x = dataFile.getVar ("wwall");
      NcVar Rwall_x = dataFile.getVar ("Rwall");
      NcVar Zwall_x = dataFile.getVar ("Zwall");
 
      wwall = new double[Nw+1];
      Rwall = new double[Nw+1];
      Zwall = new double[Nw+1];

      wwall_x.getVar (wwall);
      Rwall_x.getVar (Rwall);
      Zwall_x.getVar (Zwall);

      NcVar Rcoil_x = dataFile.getVar ("Rcoil");
      NcVar Zcoil_x = dataFile.getVar ("Zcoil");
      NcVar Icoil_x = dataFile.getVar ("Icoil");
      NcDim c_d     = Rcoil_x.getDim (0);
 
      ncoil = c_d.getSize ();

      Rcoil = new double[ncoil];
      Zcoil = new double[ncoil];
      Icoil = new double[ncoil];

      Rcoil_x.getVar (Rcoil);
      Zcoil_x.getVar (Zcoil);
      Icoil_x.getVar (Icoil);

      if (VIZ)
	{
	  NcVar RR_x = dataFile.getVar ("R");
	  NcVar ZZ_x = dataFile.getVar ("Z");
	  NcVar rr_x = dataFile.getVar ("rr");
	  NcVar tt_x = dataFile.getVar ("theta");
	  NcDim f_d  = RR_x.getDim (0);
	  
	  Nf = f_d.getSize ();
	  
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

	  delete[] RRdata; delete[] ZZdata;  delete[] rrdata; delete[] ttdata;
	}
      else
	Nf = 0;

      delete[] para;
      delete[] Hndata; delete[] Hnpdata; delete[] Vndata; delete[] Vnpdata;
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
  double* Qvac_r  = new double[J*J];
  double* Qvac_i  = new double[J*J];
  double* Svac_r  = new double[J*J];
  double* Svac_i  = new double[J*J];

  double* PRmat_r  = new double[J*J];
  double* PRmat_i  = new double[J*J];
  double* PRant_r  = new double[J*J];
  double* PRant_i  = new double[J*J];
  double* QSmat_r  = new double[J*J];
  double* QSmat_i  = new double[J*J];
  double* QSant_r  = new double[J*J];
  double* QSant_i  = new double[J*J];
  double* PSmat_r  = new double[J*J];
  double* PSmat_i  = new double[J*J];

  double* QPmat_r  = new double[J*J];
  double* QPmat_i  = new double[J*J];
  double* QPant_r  = new double[J*J];
  double* QPant_i  = new double[J*J];
  double* RSmat_r  = new double[J*J];
  double* RSmat_i  = new double[J*J];
  double* RSant_r  = new double[J*J];
  double* RSant_i  = new double[J*J];
  double* SPmat_r  = new double[J*J];
  double* SPmat_i  = new double[J*J];

  double* Hmat_r  = new double[J*J];
  double* Hmat_i  = new double[J*J];
  double* iHmat_r = new double[J*J];
  double* iHmat_i = new double[J*J];

  double* Rwal_r  = new double[J*J];
  double* Rwal_i  = new double[J*J];
  double* Swal_r  = new double[J*J];
  double* Swal_i  = new double[J*J];
  double* Imat_r  = new double[J*J];
  double* Imat_i  = new double[J*J];
  double* Iant_r  = new double[J*J];
  double* Iant_i  = new double[J*J];

  double* Gmat_r  = new double[J*J];
  double* Gmat_i  = new double[J*J];
  double* iGmat_r = new double[J*J];
  double* iGmat_i = new double[J*J];

  double* Bmat_r  = new double[J*J];
  double* Bmat_i  = new double[J*J];
  double* Cmat_r  = new double[J*J];
  double* Cmat_i  = new double[J*J];
  double* Cant_r  = new double[J*J];
  double* Cant_i  = new double[J*J];

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
  double* xii_r   = new double[J*J*NDIAG];
  double* xii_i   = new double[J*J*NDIAG];
  double* Chii_r  = new double[J*J*NDIAG];
  double* Chii_i  = new double[J*J*NDIAG];
  double* Umat_r  = new double[J*J];
  double* Umat_i  = new double[J*J];
  double* Uant_r  = new double[J*J];
  double* Uant_i  = new double[J*J];
  double* Uvec_r  = new double[J*J];
  double* Uvec_i  = new double[J*J];
  double* Psie_r  = new double[J*J*NDIAG];
  double* Psie_i  = new double[J*J*NDIAG];
  double* Ze_r    = new double[J*J*NDIAG];
  double* Ze_i    = new double[J*J*NDIAG];
  double* Xie_r   = new double[J*J*NDIAG];
  double* Xie_i   = new double[J*J*NDIAG];
  double* xie_r   = new double[J*J*NDIAG];
  double* xie_i   = new double[J*J*NDIAG];
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

  double* Wpw_r  = new double[J*J];
  double* Wpw_i  = new double[J*J];
  double* Wnw_r  = new double[J*J];
  double* Wnw_i  = new double[J*J];
  double* Wpw2_r = new double[J*J];
  double* Wpw2_i = new double[J*J];
  double* Dmat_r = new double[J*J];
  double* Dmat_i = new double[J*J];
  double* Ehmt_r = new double[J*J];
  double* Ehmt_i = new double[J*J];
  double* Fmtr_r = new double[J*J];
  double* Fmtr_i = new double[J*J];
  double* Fmta_r = new double[J*J];
  double* Fmta_i = new double[J*J];
  double* Psir_r = new double[J*J];
  double* Psir_i = new double[J*J];
  double* Xir_r  = new double[J*J];
  double* Xir_i  = new double[J*J];
  
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
	Qvac_r[cnt] = Qvac(j, jp).real();
	Qvac_i[cnt] = Qvac(j, jp).imag();
	Svac_r[cnt] = Svac(j, jp).real();
	Svac_i[cnt] = Svac(j, jp).imag();

	PRmat_r[cnt] = PRmat(j, jp).real();
	PRmat_i[cnt] = PRmat(j, jp).imag();
	PRant_r[cnt] = PRant(j, jp).real();
	PRant_i[cnt] = PRant(j, jp).imag();
	QSmat_r[cnt] = QSmat(j, jp).real();
	QSmat_i[cnt] = QSmat(j, jp).imag();
	QSant_r[cnt] = QSant(j, jp).real();
	QSant_i[cnt] = QSant(j, jp).imag();
	PSmat_r[cnt] = PSmat(j, jp).real();
	PSmat_i[cnt] = PSmat(j, jp).imag();

	QPmat_r[cnt] = QPmat(j, jp).real();
	QPmat_i[cnt] = QPmat(j, jp).imag();
	QPant_r[cnt] = QPant(j, jp).real();
	QPant_i[cnt] = QPant(j, jp).imag();
	RSmat_r[cnt] = RSmat(j, jp).real();
	RSmat_i[cnt] = RSmat(j, jp).imag();
	RSant_r[cnt] = RSant(j, jp).real();
	RSant_i[cnt] = RSant(j, jp).imag();
	SPmat_r[cnt] = SPmat(j, jp).real();
	SPmat_i[cnt] = SPmat(j, jp).imag();

	Hmat_r[cnt]  = Hmat (j, jp).real();
	Hmat_i[cnt]  = Hmat (j, jp).imag();
	iHmat_r[cnt] = iHmat(j, jp).real();
	iHmat_i[cnt] = iHmat(j, jp).imag();

	Rwal_r[cnt]  = Rwal(j, jp).real();
	Rwal_i[cnt]  = Rwal(j, jp).imag();
	Swal_r[cnt]  = Swal(j, jp).real();
	Swal_i[cnt]  = Swal(j, jp).imag();

	Imat_r[cnt]  = iImat(j, jp).real();
	Imat_i[cnt]  = iImat(j, jp).imag();
	Iant_r[cnt]  = iIant(j, jp).real();
	Iant_i[cnt]  = iIant(j, jp).imag();
	Gmat_r[cnt]  = Gmat (j, jp).real();
	Gmat_i[cnt]  = Gmat (j, jp).imag();
	iGmat_r[cnt] = iGmat(j, jp).real();
	iGmat_i[cnt] = iGmat(j, jp).imag();

	Bmat_r[cnt]  = Bmat (j, jp).real();
	Bmat_i[cnt]  = Bmat (j, jp).imag();
	Cmat_r[cnt]  = Cmat (j, jp).real();
	Cmat_i[cnt]  = Cmat (j, jp).imag();
	Cant_r[cnt]  = Cant (j, jp).real();
	Cant_i[cnt]  = Cant (j, jp).imag();
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

  if (VIZ)
    {
      int cnt = 0;
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

	if (FVAL)
	  {
	    Fvec_r[cnt] = real (Fvec(j, jp));
	    Fvec_i[cnt] = imag (Fvec(j, jp));
	  }
	cnt++;
      }
   
  for (int j = 0; j < nres; j++)
    Delta_r[j] = real (Emat(j, j));

  if (RMP)
    {
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
 
      int cnt = 0;
      for (int j = 0; j < J; j++)
	for (int i = 0; i < NDIAG; i++)
	  {
	    PPR_r[cnt] = real (Psirmp(j, i));
	    PPR_i[cnt] = imag (Psirmp(j, i));
	    ZZR_r[cnt] = real (Zrmp  (j, i));
	    ZZR_i[cnt] = imag (Zrmp  (j, i));
	    cnt++;
	  }

      if (VIZ)
	{
	  int cnt = 0;
	  for (int i = 0; i < Nf; i++)
	    for (int l = 0; l <= Nw; l++)
	      {
		PPRV_r[cnt] = real (Psirv(i, l));
		PPRV_i[cnt] = imag (Psirv(i, l));
		ZZRV_r[cnt] = real (Zrv  (i, l));
		ZZRV_i[cnt] = imag (Zrv  (i, l));
		cnt++;
	      }
	}
      
      for (int i = 0; i <= Nw; i++)
	{
	  Psis_r[i] = real (Psirmps[i]);
	  Psis_i[i] = imag (Psirmps[i]);
	  Psiz_r[i] = real (Psixs  [i]);
	  Psiz_i[i] = imag (Psixs  [i]);
	}
    }

  if (IDEAL)
    {
      int cnt = 0;
      for (int j = 0; j < J; j++)
	for (int jp = 0; jp < J; jp++)
	  for (int i = 0; i < NDIAG; i++)
	    {
	      Psii_r[cnt] = real (Psii(j, jp, i));
	      Psii_i[cnt] = imag (Psii(j, jp, i));
	      Zi_r  [cnt] = real (Zi  (j, jp, i));
	      Zi_i  [cnt] = imag (Zi  (j, jp, i));
	      Xii_r [cnt] = real (Xii (j, jp, i));
	      Xii_i [cnt] = imag (Xii (j, jp, i));
	      xii_r [cnt] = real (xii (j, jp, i));
	      xii_i [cnt] = imag (xii (j, jp, i));
	      Chii_r[cnt] = real (Chii(j, jp, i));
	      Chii_i[cnt] = imag (Chii(j, jp, i));
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
	      xie_r [cnt] = real (xie (j, jp, i));
	      xie_i [cnt] = imag (xie (j, jp, i));
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
        
      for (int j = 0; j < J; j++)
	{
	  gammax_r[j] = real (gammax[j]);
	  gammax_i[j] = imag (gammax[j]);
	  gamma_r [j] = real (gamma [j]);
	  gamma_i [j] = imag (gamma [j]);
	}

      if (RWM)
	{
	  int cnt = 0;
	  for (int j = 0; j < J; j++)
	    for (int jp = 0; jp < J; jp++)
		{
		  Wpw_r [cnt] = real (Wpw (j, jp));
		  Wpw_i [cnt] = imag (Wpw (j, jp));
		  Wnw_r [cnt] = real (Wnw (j, jp));
		  Wnw_i [cnt] = imag (Wnw (j, jp));
		  Wpw2_r[cnt] = real (Wpw2(j, jp));
		  Wpw2_i[cnt] = imag (Wpw2(j, jp));
		  Dmat_r[cnt] = real (Dmat(j, jp));
		  Dmat_i[cnt] = imag (Dmat(j, jp));
		  Ehmt_r[cnt] = real (Ehmt(j, jp));
		  Ehmt_i[cnt] = imag (Ehmt(j, jp));
		  Fmtr_r[cnt] = real (FFmt(j, jp));
		  Fmtr_i[cnt] = imag (FFmt(j, jp));
		  Fmta_r[cnt] = real (FFan(j, jp));
		  Fmta_i[cnt] = imag (FFan(j, jp));
		  Psir_r[cnt] = real (Psir(j, jp));
		  Psir_i[cnt] = imag (Psir(j, jp));
		  Xir_r [cnt] = real (Xir (j, jp));
		  Xir_i [cnt] = imag (Xir (j, jp));
		  cnt++;
		}
	}
    }
  
  try
    {
      NcFile dataFile ("../Outputs/TJ/TJ.nc", NcFile::replace);
      
      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

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
      NcVar qvacr_x = dataFile.addVar ("Qvac_r", ncDouble, vacuum_d);
      qvacr_x.putVar (Qvac_r);
      NcVar qvaci_x = dataFile.addVar ("Qvac_i", ncDouble, vacuum_d);
      qvaci_x.putVar (Qvac_i);
      NcVar svacr_x = dataFile.addVar ("Svac_r", ncDouble, vacuum_d);
      svacr_x.putVar (Svac_r);
      NcVar svaci_x = dataFile.addVar ("Svac_i", ncDouble, vacuum_d);
      svaci_x.putVar (Svac_i);
      
      NcVar prmatr_x = dataFile.addVar ("PRmat_r", ncDouble, vacuum_d);
      prmatr_x.putVar (PRmat_r);
      NcVar prmati_x = dataFile.addVar ("PRmat_i", ncDouble, vacuum_d);
      prmati_x.putVar (PRmat_i);
      NcVar prantr_x = dataFile.addVar ("PRant_r", ncDouble, vacuum_d);
      prantr_x.putVar (PRant_r);
      NcVar pranti_x = dataFile.addVar ("PRant_i", ncDouble, vacuum_d);
      pranti_x.putVar (PRant_i);
      NcVar qsvacr_x = dataFile.addVar ("QSmat_r", ncDouble, vacuum_d);
      qsvacr_x.putVar (QSmat_r);
      NcVar qsvaci_x = dataFile.addVar ("QSmat_i", ncDouble, vacuum_d);
      qsvaci_x.putVar (QSmat_i);
      NcVar qsantr_x = dataFile.addVar ("QSant_r", ncDouble, vacuum_d);
      qsantr_x.putVar (QSant_r);
      NcVar qsanti_x = dataFile.addVar ("QSant_i", ncDouble, vacuum_d);
      qsanti_x.putVar (QSant_i);
      NcVar psmatr_x = dataFile.addVar ("PSmat_r", ncDouble, vacuum_d);
      psmatr_x.putVar (PSmat_r);
      NcVar psmati_x = dataFile.addVar ("PSmat_i", ncDouble, vacuum_d);
      psmati_x.putVar (PSmat_i);

      NcVar qpmatr_x = dataFile.addVar ("QPmat_r", ncDouble, vacuum_d);
      qpmatr_x.putVar (QPmat_r);
      NcVar qpmati_x = dataFile.addVar ("QPmat_i", ncDouble, vacuum_d);
      qpmati_x.putVar (QPmat_i);
      NcVar qpantr_x = dataFile.addVar ("QPant_r", ncDouble, vacuum_d);
      qpantr_x.putVar (QPant_r);
      NcVar qpanti_x = dataFile.addVar ("QPant_i", ncDouble, vacuum_d);
      qpanti_x.putVar (QPant_i);
      NcVar rsvacr_x = dataFile.addVar ("RSmat_r", ncDouble, vacuum_d);
      rsvacr_x.putVar (RSmat_r);
      NcVar rsvaci_x = dataFile.addVar ("RSmat_i", ncDouble, vacuum_d);
      rsvaci_x.putVar (RSmat_i);
      NcVar rsantr_x = dataFile.addVar ("RSant_r", ncDouble, vacuum_d);
      rsantr_x.putVar (RSant_r);
      NcVar rsanti_x = dataFile.addVar ("RSant_i", ncDouble, vacuum_d);
      rsanti_x.putVar (RSant_i);
      NcVar spmatr_x = dataFile.addVar ("SPmat_r", ncDouble, vacuum_d);
      spmatr_x.putVar (SPmat_r);
      NcVar spmati_x = dataFile.addVar ("SPmat_i", ncDouble, vacuum_d);
      spmati_x.putVar (SPmat_i);

      NcVar hmatr_x  = dataFile.addVar ("Hmat_r",  ncDouble, vacuum_d);
      hmatr_x.putVar (Hmat_r);
      NcVar hmati_x  = dataFile.addVar ("Hmat_i",  ncDouble, vacuum_d);
      hmati_x.putVar (Hmat_i);
      NcVar ihmatr_x = dataFile.addVar ("iHmat_r", ncDouble, vacuum_d);
      ihmatr_x.putVar (iHmat_r);
      NcVar ihmati_x = dataFile.addVar ("iHmat_i", ncDouble, vacuum_d);
      ihmati_x.putVar (Hmat_i);

      NcVar rwalr_x = dataFile.addVar ("Rwal_r", ncDouble, vacuum_d);
      rwalr_x.putVar (Rwal_r);
      NcVar rwali_x = dataFile.addVar ("Rwal_i", ncDouble, vacuum_d);
      rwali_x.putVar (Rwal_i);
      NcVar swalr_x = dataFile.addVar ("Swal_r", ncDouble, vacuum_d);
      swalr_x.putVar (Swal_r);
      NcVar swali_x = dataFile.addVar ("Swal_i", ncDouble, vacuum_d);
      swali_x.putVar (Swal_i);

      NcVar imatr_x = dataFile.addVar ("iImat_r", ncDouble, vacuum_d);
      imatr_x.putVar (Imat_r);
      NcVar imati_x = dataFile.addVar ("iImat_i", ncDouble, vacuum_d);
      imati_x.putVar (Imat_i);
      NcVar iantr_x = dataFile.addVar ("iIant_r", ncDouble, vacuum_d);
      iantr_x.putVar (Iant_r);
      NcVar ianti_x = dataFile.addVar ("iIant_i", ncDouble, vacuum_d);
      ianti_x.putVar (Iant_i);

      NcVar gmatr_x  = dataFile.addVar ("Gmat_r",  ncDouble, vacuum_d);
      gmatr_x.putVar (Gmat_r);
      NcVar gmati_x  = dataFile.addVar ("Gmat_i",  ncDouble, vacuum_d);
      gmati_x.putVar (Gmat_i);
      NcVar igmatr_x = dataFile.addVar ("iGmat_r", ncDouble, vacuum_d);
      igmatr_x.putVar (iGmat_r);
      NcVar igmati_x = dataFile.addVar ("iGmat_i", ncDouble, vacuum_d);
      igmati_x.putVar (iGmat_i);

      NcVar bmatr_x  = dataFile.addVar ("Bmat_r",  ncDouble, vacuum_d);
      bmatr_x.putVar (Bmat_r);
      NcVar bmati_x  = dataFile.addVar ("Bmat_i",  ncDouble, vacuum_d);
      bmati_x.putVar (Bmat_i);
      NcVar cmatr_x  = dataFile.addVar ("Cmat_r",  ncDouble, vacuum_d);
      cmatr_x.putVar (Cmat_r);
      NcVar cmati_x  = dataFile.addVar ("Cmat_i",  ncDouble, vacuum_d);
      cmati_x.putVar (Cmat_i);
      NcVar cantr_x  = dataFile.addVar ("Cant_r",  ncDouble, vacuum_d);
      cantr_x.putVar (Cant_r);
      NcVar canti_x  = dataFile.addVar ("Cant_i",  ncDouble, vacuum_d);
      canti_x.putVar (Cant_i);
    
      NcVar rgrid_x = dataFile.addVar ("r_grid",    ncDouble, d_d);
      rgrid_x.putVar (Rgrid);
      NcVar pgrid_x = dataFile.addVar ("PsiN_grid", ncDouble, d_d);
      pgrid_x.putVar (Pgrid);
      NcVar hode_x  = dataFile.addVar ("h_ode",     ncDouble, d_d);
      hode_x.putVar (hode);
      NcVar eode_x  = dataFile.addVar ("err_ode",   ncDouble, d_d);
      eode_x.putVar (eode);

      NcVar t_x      = dataFile.addVar ("theta",  ncDouble, w_d);
      t_x.putVar (tbound);
      NcVar cmu_x    = dataFile.addVar ("cosmu",  ncDouble, w_d);
      cmu_x.putVar (cmu);
      NcVar e_x      = dataFile.addVar ("eta",    ncDouble, w_d);
      e_x.putVar (eeta);
      NcVar ceta_x   = dataFile.addVar ("coseta", ncDouble, w_d);
      ceta_x.putVar (ceta);
      NcVar seta_x   = dataFile.addVar ("sineta", ncDouble, w_d);
      seta_x.putVar (seta);
      NcVar R2grgz_x = dataFile.addVar ("R2grgz", ncDouble, w_d);
      R2grgz_x.putVar (R2grgz);
      NcVar R2grge_x = dataFile.addVar ("R2grge", ncDouble, w_d);
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

      NcVar mres_x   = dataFile.addVar ("m_res",    ncInt,    x_d);
      mres_x.putVar (mres);
      NcVar rres_x   = dataFile.addVar ("r_res",    ncDouble, x_d);
      rres_x.putVar (rres);
      NcVar sres_x   = dataFile.addVar ("s_res",    ncDouble, x_d);
      sres_x.putVar (sres);
      NcVar dires_x  = dataFile.addVar ("DI_res",   ncDouble, x_d);
      dires_x.putVar (DIres);
      NcVar drres_x  = dataFile.addVar ("DR_res",   ncDouble, x_d);
      drres_x.putVar (DRres);
      NcVar flarge_x = dataFile.addVar ("Flarge",   ncDouble, x_d);
      flarge_x.putVar (Flarge);
      NcVar fsmall_x = dataFile.addVar ("Fsmall",   ncDouble, x_d);
      fsmall_x.putVar (Fsmall);
      NcVar Pn_x     = dataFile.addVar ("PsiN_res", ncDouble, x_d);
      Pn_x.putVar (Pres);

      NcVar tf_x = dataFile.addVar ("Torque_full", ncDouble, t_d);
      tf_x.putVar (Tf.data());
      NcVar tu_x = dataFile.addVar ("Torque_unrc", ncDouble, t_d);
      tu_x.putVar (Tu.data());

      NcVar tfull_x = dataFile.addVar ("Torque_pair_full", ncDouble, tt_d);
      tfull_x.putVar (Tfull.data());
      NcVar tunrc_x = dataFile.addVar ("Torque_pair_unrc", ncDouble, tt_d);
      tunrc_x.putVar (Tunrc.data());

      if (VIZ)
	{
	  NcVar ppvr_x = dataFile.addVar ("Psi_unrc_eig_r", ncDouble, v_d);
	  ppvr_x.putVar (PPV_r);
	  NcVar ppvi_x = dataFile.addVar ("Psi_unrc_eig_i", ncDouble, v_d);
	  ppvi_x.putVar (PPV_i);
	  NcVar zzvr_x = dataFile.addVar ("Z_unrc_eig_r",   ncDouble, v_d);
	  zzvr_x.putVar (ZZV_r);
	  NcVar zzvi_x = dataFile.addVar ("Z_unrc_eig_i",   ncDouble, v_d);
	  zzvi_x.putVar (ZZV_i);
	}

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

      if (FVAL)
	{
	  NcVar Fval_x  = dataFile.addVar ("Fval",   ncDouble, x_d);
	  Fval_x.putVar (Fval);
	  NcVar Fvecr_x = dataFile.addVar ("Fvec_r", ncDouble, e_d);
	  Fvecr_x.putVar (Fvec_r);
	  NcVar Fveci_x = dataFile.addVar ("Fvec_i", ncDouble, e_d);
	  Fveci_x.putVar (Fvec_i);
	}

      NcVar rcoil_x = dataFile.addVar ("Rcoil",     ncDouble, c_d);
      rcoil_x.putVar (Rcoil);
      NcVar zcoil_x = dataFile.addVar ("Zcoil",     ncDouble, c_d);
      zcoil_x.putVar (Zcoil);
      NcVar icoil_x = dataFile.addVar ("Icoil",     ncDouble, c_d);
      icoil_x.putVar (Icoil);

      if (RMP)
	{
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

	  NcVar pprr_x = dataFile.addVar ("Psi_rmp_r", ncDouble, rmp_d);
	  pprr_x.putVar (PPR_r);
	  NcVar ppri_x = dataFile.addVar ("Psi_rmp_i", ncDouble, rmp_d);
	  ppri_x.putVar (PPR_i);
	  NcVar zzrr_x = dataFile.addVar ("Z_rmp_r",   ncDouble, rmp_d);
	  zzrr_x.putVar (ZZR_r);
	  NcVar zzri_x = dataFile.addVar ("Z_rmp_i",   ncDouble, rmp_d);
	  zzri_x.putVar (ZZR_i);

	  if (VIZ)
	    {
	      NcVar pprvr_x = dataFile.addVar ("Psi_rmp_eig_r", ncDouble, vr_d);
	      pprvr_x.putVar (PPRV_r);
	      NcVar pprvi_x = dataFile.addVar ("Psi_rmp_eig_i", ncDouble, vr_d);
	      pprvi_x.putVar (PPRV_i);
	      NcVar zzrvr_x = dataFile.addVar ("Z_rmp_eig_r",   ncDouble, vr_d);
	      zzrvr_x.putVar (ZZRV_r);
	      NcVar zzrvi_x = dataFile.addVar ("Z_rmp_eig_i",   ncDouble, vr_d);
	      zzrvi_x.putVar (ZZRV_i);
	    }
	  
	  NcVar psisr_x = dataFile.addVar ("Psi_rmp_surface_r", ncDouble, w_d);
	  psisr_x.putVar (Psis_r);
	  NcVar psisi_x = dataFile.addVar ("Psi_rmp_surface_i", ncDouble, w_d);
	  psisi_x.putVar (Psis_i);
	  NcVar psixr_x = dataFile.addVar ("Psi_x_surface_r",   ncDouble, w_d);
	  psixr_x.putVar (Psiz_r);
	  NcVar psixi_x = dataFile.addVar ("Psi_x_surface_i",   ncDouble, w_d);
	  psixi_x.putVar (Psiz_i);
	}

      if (IDEAL)
	{
	  NcVar psiir_x = dataFile.addVar ("Psi_i_r", ncDouble, ideal_d);
	  psiir_x.putVar (Psii_r);
	  NcVar psiii_x = dataFile.addVar ("Psi_i_i", ncDouble, ideal_d);
	  psiii_x.putVar (Psii_i);
	  NcVar zir_x   = dataFile.addVar ("Z_i_r",   ncDouble, ideal_d);
	  zir_x.putVar (Zi_r);
	  NcVar zii_x   = dataFile.addVar ("Z_i_i",   ncDouble, ideal_d);
	  zii_x.putVar (Zi_i);
	  NcVar Xir_x   = dataFile.addVar ("Xi_i_r",  ncDouble, ideal_d);
	  Xir_x.putVar (Xii_r);
	  NcVar Xii_x   = dataFile.addVar ("Xi_i_i",  ncDouble, ideal_d);
	  Xii_x.putVar (Xii_i);
	  NcVar xir_x   = dataFile.addVar ("xi_i_r",  ncDouble, ideal_d);
	  xir_x.putVar (xii_r);
	  NcVar xii_x   = dataFile.addVar ("xi_i_i",  ncDouble, ideal_d);
	  xii_x.putVar (xii_i);
	  NcVar chiir_x = dataFile.addVar ("Chi_i_r", ncDouble, ideal_d);
	  chiir_x.putVar (Chii_r);
	  NcVar chiii_x = dataFile.addVar ("Chi_i_i", ncDouble, ideal_d);
	  chiii_x.putVar (Chii_i);
	  if (INTR)
	    {
	      NcVar lvals_x = dataFile.addVar ("lambda",  ncDouble, rmp_d);
	      lvals_x.putVar (xvals);
	    }
	  
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
	  
	  NcVar psier_x = dataFile.addVar ("Psi_e_r", ncDouble, ideal_d);
	  psier_x.putVar (Psie_r);
	  NcVar psiei_x = dataFile.addVar ("Psi_e_i", ncDouble, ideal_d);
	  psiei_x.putVar (Psie_i);
	  NcVar zer_x   = dataFile.addVar ("Z_e_r",   ncDouble, ideal_d);
	  zer_x.putVar (Ze_r);
	  NcVar zei_x   = dataFile.addVar ("Z_e_i",   ncDouble, ideal_d);
	  zei_x.putVar (Ze_i);
	  NcVar Xer_x   = dataFile.addVar ("Xi_e_r",  ncDouble, ideal_d);
	  Xer_x.putVar (Xie_r);
	  NcVar Xei_x   = dataFile.addVar ("Xi_e_i",  ncDouble, ideal_d);
	  Xei_x.putVar (Xie_i);
	  NcVar xer_x   = dataFile.addVar ("xi_e_r",  ncDouble, ideal_d);
	  xer_x.putVar (xie_r);
	  NcVar xei_x   = dataFile.addVar ("xi_e_i",  ncDouble, ideal_d);
	  xei_x.putVar (xie_i);
	  
	  
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

	  if (RMP)
	    {
	      NcVar gammaxr_x = dataFile.addVar ("gammax_r", ncDouble, j_d);
	      gammaxr_x.putVar (gammax_r);
	      NcVar gammaxi_x = dataFile.addVar ("gammax_i", ncDouble, j_d);
	      gammaxi_x.putVar (gammax_i);
	      NcVar gammar_x = dataFile.addVar  ("gamma_r",  ncDouble, j_d);
	      gammar_x.putVar (gamma_r);
	      NcVar gammai_x = dataFile.addVar  ("gamma_i",  ncDouble, j_d);
	      gammai_x.putVar (gamma_i);
	    }

	  if (RWM)
	    {
	      NcVar wpwr_x  = dataFile.addVar ("Wpw_r",     ncDouble, vacuum_d);
	      wpwr_x.putVar (Wpw_r);
	      NcVar wpwi_x  = dataFile.addVar ("Wpw_i",     ncDouble, vacuum_d);
	      wpwi_x.putVar (Wpw_i);
	      NcVar wnwr_x  = dataFile.addVar ("Wnw_r",     ncDouble, vacuum_d);
	      wnwr_x.putVar (Wnw_r);
	      NcVar wnwi_x  = dataFile.addVar ("Wnw_i",     ncDouble, vacuum_d);
	      wnwi_x.putVar (Wnw_i);
	      NcVar wpw2r_x = dataFile.addVar ("Wpw2_r",    ncDouble, vacuum_d);
	      wpw2r_x.putVar (Wpw2_r);
	      NcVar wpw2i_x = dataFile.addVar ("Wpw2_i",    ncDouble, vacuum_d);
	      wpw2i_x.putVar (Wpw2_i);
	      NcVar dmatr_x = dataFile.addVar ("Dmat_r",    ncDouble, vacuum_d);
	      dmatr_x.putVar (Dmat_r);
	      NcVar dmati_x = dataFile.addVar ("Dmat_i",    ncDouble, vacuum_d);
	      dmati_x.putVar (Dmat_i);
	      NcVar ehmtr_x = dataFile.addVar ("Ehmt_r",    ncDouble, vacuum_d);
	      ehmtr_x.putVar (Ehmt_r);
	      NcVar ehmti_x = dataFile.addVar ("Ehmt_i",    ncDouble, vacuum_d);
	      ehmti_x.putVar (Ehmt_i);
	      NcVar fmtrr_x = dataFile.addVar ("Fmtr_r",    ncDouble, vacuum_d);
	      fmtrr_x.putVar (Fmtr_r);
	      NcVar fmtri_x = dataFile.addVar ("Fmtr_i",    ncDouble, vacuum_d);
	      fmtri_x.putVar (Fmtr_i);
	      NcVar fmtar_x = dataFile.addVar ("Fmta_r",    ncDouble, vacuum_d);
	      fmtar_x.putVar (Fmta_r);
	      NcVar fmtai_x = dataFile.addVar ("Fmta_i",    ncDouble, vacuum_d);
	      fmtai_x.putVar (Fmta_i);
	      NcVar fw_x    = dataFile.addVar ("f_w",       ncDouble, j_d);
	      fw_x.putVar (fw);
	      NcVar psirr_x = dataFile.addVar ("Psirwm_r",  ncDouble, vacuum_d);
	      psirr_x.putVar (Psir_r);
	      NcVar psiri_x = dataFile.addVar ("Psirwm_i",  ncDouble, vacuum_d);
	      psiri_x.putVar (Psir_i);
	      NcVar xirr_x  = dataFile.addVar ("Xi_rwm_r",  ncDouble, vacuum_d);
	      xirr_x.putVar (Xir_r);
	      NcVar xiri_x  = dataFile.addVar ("Xi_rwm_i",  ncDouble, vacuum_d);
	      xiri_x.putVar (Xir_i);
	    }
	}

      if (LAYER)
	{
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
	}  
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

  delete[] Pvac_r;  delete[] Pvac_i;  delete[] Rvac_r;  delete[] Rvac_i;
  delete[] PRmat_r; delete[] PRmat_i; delete[] PRant_r; delete[] PRant_i;
  delete[] QSmat_r; delete[] QSmat_i; delete[] QSant_r; delete[] QSant_i;
  delete[] Hmat_r;  delete[] Hmat_i;  delete[] Qvac_r;  delete[] Qvac_i;
  delete[] Svac_r;  delete[] Svac_i;  delete[] PSmat_r; delete[] PSmat_i;
  delete[] iHmat_r; delete[] iHmat_i; delete[] SPmat_r; delete[] SPmat_i;
  delete[] QPmat_r; delete[] QPmat_i; delete[] QPant_r; delete[] QPant_i;
  delete[] RSmat_r; delete[] RSmat_i; delete[] RSant_r; delete[] RSant_i;
  delete[] Rwal_r;  delete[] Rwal_i;  delete[] Swal_r;  delete[] Swal_i;
  delete[] Imat_r;  delete[] Imat_i;  delete[] Gmat_r;  delete[] Gmat_i;
  delete[] iGmat_r; delete[] iGmat_i; delete[] Bmat_r;  delete[] Bmat_i;
  delete[] Cmat_r;  delete[] Cmat_i;  delete[] Cant_r;  delete[] Cant_i;
  delete[] Iant_r;  delete[] Iant_i;  
   
  delete[] PPPsi_r; delete[] PPPsi_i; delete[] ZZZ_r; delete[] ZZZ_i;
  delete[] PPF_r;   delete[] PPF_i;   delete[] ZZF_r; delete[] ZZF_i;
  delete[] PPU_r;   delete[] PPU_i;   delete[] ZZU_r; delete[] ZZU_i;

  delete[] PPV_r;   delete[] PPV_i;  delete[] ZZV_r;  delete[] ZZV_i;
  delete[] Emat_r;  delete[] Emat_i; delete[] Eant_r; delete[] Eant_i;
  delete[] Psix_r;  delete[] Psix_i; delete[] Xi_r;   delete[] Xi_i; 
  delete[] Up_r;    delete[] Up_i;   delete[] Chi_r;  delete[] Chi_i;
  delete[] Fvec_r;  delete[] Fvec_i; delete[] Fmat_r; delete[] Fmat_i;
  delete[] Delta_r;

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
  delete[] Chii_r;   delete[] Chii_i;   delete[] xvals;   delete[] xii_r;
  delete[] xii_i;    delete[] xie_r;    delete[] xie_i; 
 
  delete[] Wpw_r;  delete[] Wpw_i;  delete[] Wnw_r;  delete[] Wnw_i;
  delete[] Wpw2_r; delete[] Wpw2_i; delete[] Dmat_r; delete[] Dmat_i;
  delete[] Fmtr_r; delete[] Fmtr_i; delete[] Fmta_r; delete[] Fmta_i;
  delete[] Ehmt_r; delete[] Ehmt_i; delete[] Psir_r; delete[] Psir_i;
  delete[] Xir_r;  delete[] Xir_i;
}
