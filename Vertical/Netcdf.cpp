// Netcdf.cpp

#include "Vertical.h"

#define NINPUT 21
#define NOUTPUT 2

// ##################################################
// Function to read equilibrium data from netcdf file
// ##################################################
void Vertical::ReadEquilibriumNetcdf ()
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
      H1b  = para[3];
      H1pb = para[4];
      
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
      NcVar S4_x  = dataFile.getVar ("S4");
      NcVar S5_x  = dataFile.getVar ("S5");
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P3_x  = dataFile.getVar ("P3");
      NcVar P4_x  = dataFile.getVar ("P3a");
      NcVar ne_x  = dataFile.getVar ("ne");
      NcVar Te_x  = dataFile.getVar ("Te");
      NcVar nep_x = dataFile.getVar ("nep");
      NcVar Tep_x = dataFile.getVar ("Tep");
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
      S4   = new double[Nr+1];
      S5   = new double[Nr+1];
      P1   = new double[Nr+1];
      P2   = new double[Nr+1];
      P3   = new double[Nr+1];
      Sig  = new double[Nr+1];
      ne   = new double[Nr+1];
      Te   = new double[Nr+1];
      nep  = new double[Nr+1];
      Tep  = new double[Nr+1];

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
      S4_x .getVar (S4);
      S5_x .getVar (S5);
      P1_x .getVar (P1);
      P2_x .getVar (P2);
      P3_x .getVar (P3);
      P4_x .getVar (Sig);
      ne_x .getVar (ne);
      Te_x .getVar (Te);
      nep_x.getVar (nep);
      Tep_x.getVar (Tep);

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

      NcVar twall_x  = dataFile.getVar ("twall");
      NcVar Rwall_x  = dataFile.getVar ("Rwall");
      NcVar Zwall_x  = dataFile.getVar ("Zwall");
      NcVar dRdthw_x = dataFile.getVar ("dRdthetawall");
      NcVar dZdthw_x = dataFile.getVar ("dZdthetawall");
 
      twall  = new double[Nw+1];
      Rwall  = new double[Nw+1];
      Zwall  = new double[Nw+1];
      dRdthw = new double[Nw+1];
      dZdthw = new double[Nw+1];

      twall_x .getVar (twall);
      Rwall_x .getVar (Rwall);
      Zwall_x .getVar (Zwall);
      dRdthw_x.getVar (dRdthw);
      dZdthw_x.getVar (dZdthw);

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

	  delete[] RRdata; delete[] ZZdata; delete[] rrdata; delete[] ttdata;

	  dataFile.close ();
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
void Vertical::WriteNetcdf ()
{
  printf ("Writing stability data to netcdf file Outputs/Vertical/Vertical.nc:\n");

  double Input [NINPUT];
  double Output[NOUTPUT];

  double* Pvac_r = new double[J*J];
  double* Pvac_i = new double[J*J];
  double* Rvac_r = new double[J*J];
  double* Rvac_i = new double[J*J];
  double* Qvac_r = new double[J*J];
  double* Qvac_i = new double[J*J];
  double* Svac_r = new double[J*J];
  double* Svac_i = new double[J*J];

  double* PRmat_r = new double[J*J];
  double* PRmat_i = new double[J*J];
  double* PRant_r = new double[J*J];
  double* PRant_i = new double[J*J];
  double* QSmat_r = new double[J*J];
  double* QSmat_i = new double[J*J];
  double* QSant_r = new double[J*J];
  double* QSant_i = new double[J*J];
  double* PSmat_r = new double[J*J];
  double* PSmat_i = new double[J*J];

  double* QPmat_r = new double[J*J];
  double* QPmat_i = new double[J*J];
  double* QPant_r = new double[J*J];
  double* QPant_i = new double[J*J];
  double* RSmat_r = new double[J*J];
  double* RSmat_i = new double[J*J];
  double* RSant_r = new double[J*J];
  double* RSant_i = new double[J*J];
  double* SPmat_r = new double[J*J];
  double* SPmat_i = new double[J*J];

  double* Hmat_r  = new double[J*J];
  double* Hmat_i  = new double[J*J];
  double* iHmat_r = new double[J*J];
  double* iHmat_i = new double[J*J];

  double* Rwal_r = new double[J*J];
  double* Rwal_i = new double[J*J];
  double* Swal_r = new double[J*J];
  double* Swal_i = new double[J*J];
  double* Imat_r = new double[J*J];
  double* Imat_i = new double[J*J];
  double* Iant_r = new double[J*J];
  double* Iant_i = new double[J*J];

  double* Gmat_r  = new double[J*J];
  double* Gmat_i  = new double[J*J];
  double* iGmat_r = new double[J*J];
  double* iGmat_i = new double[J*J];

  double* Ammp_r = new double[(Nr+1)*J*J];
  double* Bmmp_r = new double[(Nr+1)*J*J];
  double* Cmmp_r = new double[(Nr+1)*J*J];
  double* Dmmp_r = new double[(Nr+1)*J*J];
  double* Ammp_i = new double[(Nr+1)*J*J];
  double* Bmmp_i = new double[(Nr+1)*J*J];
  double* Cmmp_i = new double[(Nr+1)*J*J];
  double* Dmmp_i = new double[(Nr+1)*J*J];
  double* Atest  = new double[(Nr+1)*J*J];
  double* BCtest = new double[(Nr+1)*J*J];
  double* Dtest  = new double[(Nr+1)*J*J];

  double* PPPsi_r = new double[J*J*NDIAG];
  double* PPPsi_i = new double[J*J*NDIAG];
  double* ZZZ_r   = new double[J*J*NDIAG];
  double* ZZZ_i   = new double[J*J*NDIAG];

  double* Umat_r = new double[J*J];
  double* Umat_i = new double[J*J];
  double* Uant_r = new double[J*J];
  double* Uant_i = new double[J*J];
  double* Uvec_r = new double[J*J];
  double* Uvec_i = new double[J*J];
  double* Psie_r = new double[J*J*NDIAG];
  double* Psie_i = new double[J*J*NDIAG];
  double* Ze_r   = new double[J*J*NDIAG];
  double* Ze_i   = new double[J*J*NDIAG];
  
  double* Psiy_r = new double[J*(Nw+1)];
  double* Psiy_i = new double[J*(Nw+1)];
  double* Xiy_r  = new double[J*(Nw+1)];
  double* Xiy_i  = new double[J*(Nw+1)];

  double* pUmat_r = new double[J*J];
  double* pUmat_i = new double[J*J];
  double* pUant_r = new double[J*J];
  double* pUant_i = new double[J*J];
  double* pUvec_r = new double[J*J];
  double* pUvec_i = new double[J*J];
  double* pPsie_r = new double[J*J*NDIAG];
  double* pPsie_i = new double[J*J*NDIAG];
  double* pZe_r   = new double[J*J*NDIAG];
  double* pZe_i   = new double[J*J*NDIAG];
  
  double* pPsiy_r = new double[J*(Nw+1)];
  double* pPsiy_i = new double[J*(Nw+1)];
  double* pXiy_r  = new double[J*(Nw+1)];
  double* pXiy_i  = new double[J*(Nw+1)];

  double* PPV_r = new double[J*Nf*(Nw+1)];
  double* PPV_i = new double[J*Nf*(Nw+1)];
  double* ZZV_r = new double[J*Nf*(Nw+1)];
  double* ZZV_i = new double[J*Nf*(Nw+1)];

  double* pPPV_r = new double[J*Nf*(Nw+1)];
  double* pPPV_i = new double[J*Nf*(Nw+1)];
  double* pZZV_r = new double[J*Nf*(Nw+1)];
  double* pZZV_i = new double[J*Nf*(Nw+1)];

  double* ya_r = new double[J];
  double* ya_i = new double[J];
  double* yb_r = new double[J];
  double* yb_i = new double[J];
  
  Input[0]  = double (MMIN);
  Input[1]  = double (MMAX);
  Input[2]  = EPS;
  Input[3]  = double (NFIX);
  Input[4]  = double (NDIAG);
  Input[5]  = acc;
  Input[6]  = h0;
  Input[7]  = hmin;
  Input[8]  = hmax;
  Input[9]  = EPSF;
  Input[11] = double (SRC);
  Input[12] = B0;
  Input[13] = R0;
  Input[14] = n0;
  Input[15] = alpha;
  Input[16] = Zeff;
  Input[17] = Mion;
  Input[18] = Chip;
  Input[19] = Teped;
  Input[20] = neped;

  Output[0] = alphaw;
  Output[1] = gammaw;

  int cnt = 0;
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

	Hmat_r [cnt] = Hmat (j, jp).real();
	Hmat_i [cnt] = Hmat (j, jp).imag();
	iHmat_r[cnt] = iHmat(j, jp).real();
	iHmat_i[cnt] = iHmat(j, jp).imag();
	
	Rwal_r[cnt] = Rwal(j, jp).real();
	Rwal_i[cnt] = Rwal(j, jp).imag();
	Swal_r[cnt] = Swal(j, jp).real();
	Swal_i[cnt] = Swal(j, jp).imag();

	Imat_r [cnt] = iImat(j, jp).real();
	Imat_i [cnt] = iImat(j, jp).imag();
	Iant_r [cnt] = iIant(j, jp).real();
	Iant_i [cnt] = iIant(j, jp).imag();
	Gmat_r [cnt] = Gmat (j, jp).real();
	Gmat_i [cnt] = Gmat (j, jp).imag();
	iGmat_r[cnt] = iGmat(j, jp).real();
	iGmat_i[cnt] = iGmat(j, jp).imag();

	cnt++;
      }
  
  cnt = 0;
  for (int i = 0; i <= Nr; i++)
    for (int j = 0; j < J; j++)
      for (int jp = 0; jp < J; jp++)
	{
	  complex<double> ammp, Ammp, Bmmp, Cmmp, Dmmp;
	  complex<double> ampm, Ampm, Bmpm, Cmpm, Dmpm;
	  GetMatrices (rr[i], MPOL[j],  MPOL[jp], Ammp, Bmmp, Cmmp, Dmmp);
	  GetMatrices (rr[i], MPOL[jp], MPOL[j],  Ampm, Bmpm, Cmpm, Dmpm);
	  
	  Ammp_r[cnt] = Ammp.real();
	  Ammp_i[cnt] = Ammp.imag();
	  Bmmp_r[cnt] = Bmmp.real();
	  Bmmp_i[cnt] = Bmmp.imag();
	  Cmmp_r[cnt] = Cmmp.real();
	  Cmmp_i[cnt] = Cmmp.imag();
	  Dmmp_r[cnt] = Dmmp.real();
	  Dmmp_i[cnt] = Dmmp.imag();

	  Atest [cnt] = abs (Ammp - conj (Ampm));
	  BCtest[cnt] = abs (Bmmp + conj (Cmpm));
	  Dtest [cnt] = abs (Dmmp - conj (Dmpm));
	  
	  cnt++;
	}

  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
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
    for (int jp = 0; jp < J; jp++)
      {
	Umat_r[cnt]  = Umat(j, jp) .real();
	Umat_i[cnt]  = Umat(j, jp) .imag();
	Uant_r[cnt]  = Uant(j, jp) .real();
	Uant_i[cnt]  = Uant(j, jp) .imag();
	Uvec_r[cnt]  = Uvec(j, jp) .real();
	Uvec_i[cnt]  = Uvec(j, jp) .imag();

	pUmat_r[cnt]  = pUmat(j, jp) .real();
	pUmat_i[cnt]  = pUmat(j, jp) .imag();
	pUant_r[cnt]  = pUant(j, jp) .real();
	pUant_i[cnt]  = pUant(j, jp) .imag();
	pUvec_r[cnt]  = pUvec(j, jp) .real();
	pUvec_i[cnt]  = pUvec(j, jp) .imag();
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

	  pPsie_r[cnt] = real (pPsie(j, jp, i));
	  pPsie_i[cnt] = imag (pPsie(j, jp, i));
	  pZe_r  [cnt] = real (pZe  (j, jp, i));
	  pZe_i  [cnt] = imag (pZe  (j, jp, i));
	  cnt++;
	}
      
  cnt = 0;
  for (int j = 0; j < J; j++)
    for (int i = 0; i <= Nw; i++)
      {
	Psiy_r[cnt] = real (Psiy(j, i));
	Psiy_i[cnt] = imag (Psiy(j, i));
	Xiy_r [cnt] = real (Xiy (j, i));
	Xiy_i [cnt] = imag (Xiy (j, i));

	pPsiy_r[cnt] = real (pPsiy(j, i));
	pPsiy_i[cnt] = imag (pPsiy(j, i));
	pXiy_r [cnt] = real (pXiy (j, i));
	pXiy_i [cnt] = imag (pXiy (j, i));
	cnt++;
      }

  for (int j = 0; j < J; j++)
    {
      ya_r[j] = real (ya[j]);
      ya_i[j] = imag (ya[j]);
      yb_r[j] = real (yb[j]);
      yb_i[j] = imag (yb[j]);
    }

  if (VIZ)
    {
      int cnt = 0;
      for (int k = 0; k < J; k++)
	for (int i = 0; i < Nf; i++)
	  for (int l = 0; l <= Nw; l++)
	    {
	      PPV_r[cnt] = real (Psiuv(k, i, l));
	      PPV_i[cnt] = imag (Psiuv(k, i, l));
	      ZZV_r[cnt] = real (Zuv  (k, i, l));
	      ZZV_i[cnt] = imag (Zuv  (k, i, l));

	      pPPV_r[cnt] = real (pPsiuv(k, i, l));
	      pPPV_i[cnt] = imag (pPsiuv(k, i, l));
	      pZZV_r[cnt] = real (pZuv  (k, i, l));
	      pZZV_i[cnt] = imag (pZuv  (k, i, l));
	      cnt++;
	    }
    }
  
  try
    {
      NcFile dataFile ("../Outputs/Vertical/Vertical.nc", NcFile::replace);

      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim i_d = dataFile.addDim ("Ni",    NINPUT);
      NcDim o_d = dataFile.addDim ("No",    NOUTPUT);
      NcDim r_d = dataFile.addDim ("Nr",    Nr+1);
      NcDim s_d = dataFile.addDim ("Ns",    Ns+1);
      NcDim j_d = dataFile.addDim ("J",     J);
      NcDim w_d = dataFile.addDim ("Nw",    Nw+1);
      NcDim d_d = dataFile.addDim ("ndiag", NDIAG);
      NcDim f_d = dataFile.addDim ("Nf",    Nf);

      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> vacuum_d;
      vacuum_d.push_back (j_d);
      vacuum_d.push_back (j_d);

      vector<NcDim> matrix_d;
      matrix_d.push_back (r_d);
      matrix_d.push_back (j_d);
      matrix_d.push_back (j_d);

      vector<NcDim> soln_d;
      soln_d.push_back (j_d);
      soln_d.push_back (j_d);
      soln_d.push_back (d_d);

      vector<NcDim> torque_d;
      torque_d.push_back (j_d);
      torque_d.push_back (d_d);

      vector<NcDim> surface_d;
      surface_d.push_back (j_d);
      surface_d.push_back (w_d);
      
      vector<NcDim> ideal_d;
      ideal_d.push_back (j_d);
      ideal_d.push_back (j_d);
      ideal_d.push_back (d_d);
      
      vector<NcDim> v_d;
      v_d.push_back (j_d);
      v_d.push_back (f_d);
      v_d.push_back (w_d);

      NcVar i_x = dataFile.addVar ("InputParameters",  ncDouble, i_d);
      i_x.putVar (Input);
      NcVar o_x = dataFile.addVar ("OutputParameters", ncDouble, o_d);
      o_x.putVar (Output);

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
      NcVar S5_x  = dataFile.addVar ("S5",  ncDouble, r_d);
      S5_x.putVar (S5);
      NcVar Sig_x = dataFile.addVar ("Sig", ncDouble, r_d);
      Sig_x.putVar (Sig);
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

      NcVar cmuw_x    = dataFile.addVar ("cosmuwall",  ncDouble, w_d);
      cmuw_x.putVar (cmuw);
      NcVar ew_x      = dataFile.addVar ("etawall",    ncDouble, w_d);
      ew_x.putVar (eetaw);
      NcVar cetaw_x   = dataFile.addVar ("cosetawall", ncDouble, w_d);
      cetaw_x.putVar (cetaw);
      NcVar setaw_x   = dataFile.addVar ("sinetawall", ncDouble, w_d);
      setaw_x.putVar (setaw);
      NcVar R2grgzw_x = dataFile.addVar ("R2grgzwall", ncDouble, w_d);
      R2grgzw_x.putVar (R2grgzw);
      NcVar R2grgew_x = dataFile.addVar ("R2grgewall", ncDouble, w_d);
      R2grgew_x.putVar (R2grgew);
      
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

      NcVar ammpr_x = dataFile.addVar ("Ammp_r", ncDouble, matrix_d);
      ammpr_x.putVar (Ammp_r);
      NcVar ammpi_x = dataFile.addVar ("Ammp_i", ncDouble, matrix_d);
      ammpi_x.putVar (Ammp_i);
      
      NcVar bmmpr_x = dataFile.addVar ("Bmmp_r", ncDouble, matrix_d);
      bmmpr_x.putVar (Bmmp_r);
      NcVar bmmpi_x = dataFile.addVar ("Bmmp_i", ncDouble, matrix_d);
      bmmpi_x.putVar (Bmmp_i);
      
      NcVar cmmpr_x = dataFile.addVar ("Cmmp_r", ncDouble, matrix_d);
      cmmpr_x.putVar (Cmmp_r);
      NcVar cmmpi_x = dataFile.addVar ("Cmmp_i", ncDouble, matrix_d);
      cmmpi_x.putVar (Cmmp_i);
      
      NcVar dmmpr_x = dataFile.addVar ("Dmmp_r", ncDouble, matrix_d);
      dmmpr_x.putVar (Dmmp_r);
      NcVar dmmpi_x = dataFile.addVar ("Dmmp_i", ncDouble, matrix_d);
      dmmpi_x.putVar (Dmmp_i);

      NcVar atest_x  = dataFile.addVar ("Atest",  ncDouble, matrix_d);
      atest_x.putVar (Atest);
      NcVar bctest_x = dataFile.addVar ("BCtest", ncDouble, matrix_d);
      bctest_x.putVar (BCtest);
      NcVar dtest_x  = dataFile.addVar ("Dtest",  ncDouble, matrix_d);
      dtest_x.putVar (Dtest);

      NcVar pppsir_x = dataFile.addVar ("y_r", ncDouble, soln_d);
      pppsir_x.putVar (PPPsi_r);
      NcVar pppsii_x = dataFile.addVar ("y_i", ncDouble, soln_d);
      pppsii_x.putVar (PPPsi_i);
      NcVar zzzr_x   = dataFile.addVar ("Z_r", ncDouble, soln_d);
      zzzr_x.putVar (ZZZ_r);
      NcVar zzzi_x   = dataFile.addVar ("Z_i",  ncDouble, soln_d);
      zzzi_x.putVar (ZZZ_i);

      NcVar rgrid_x = dataFile.addVar ("r_grid",    ncDouble, d_d);
      rgrid_x.putVar (Rgrid);
      NcVar pgrid_x = dataFile.addVar ("PsiN_grid", ncDouble, d_d);
      pgrid_x.putVar (Pgrid);
      NcVar hode_x  = dataFile.addVar ("h_ode",     ncDouble, d_d);
      hode_x.putVar (hode);
      NcVar eode_x  = dataFile.addVar ("err_ode",   ncDouble, d_d);
      eode_x.putVar (eode);

      NcVar ttest_x = dataFile.addVar ("Energy_test", ncDouble, torque_d);
      ttest_x.putVar (Etest.data());
      NcVar pnorm_x = dataFile.addVar ("y_norm",      ncDouble, torque_d);
      pnorm_x.putVar (Pnorm.data());
      NcVar znorm_x = dataFile.addVar ("Z_norm",      ncDouble, torque_d);
      znorm_x.putVar (Znorm.data());

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

      NcVar pwvacr_x  = dataFile.addVar ("pUmat_r",  ncDouble, vacuum_d);
      pwvacr_x.putVar (pUmat_r);
      NcVar pwvaci_x  = dataFile.addVar ("pUmat_i",  ncDouble, vacuum_d);
      pwvaci_x.putVar (pUmat_i);
      NcVar pwantr_x  = dataFile.addVar ("pUant_r",  ncDouble, vacuum_d);
      pwantr_x.putVar (pUant_r);
      NcVar pwanti_x  = dataFile.addVar ("pUant_i",  ncDouble, vacuum_d);
      pwanti_x.putVar (pUant_i);
      NcVar pwval_x   = dataFile.addVar ("pWval",    ncDouble, j_d);
      pwval_x.putVar  (pWval);
      NcVar pvval_x   = dataFile.addVar ("pVval",    ncDouble, j_d);
      pvval_x.putVar  (pVval);
      NcVar puval_x   = dataFile.addVar ("pUval",    ncDouble, j_d);
      puval_x.putVar  (pUval);
      NcVar pwvecr_x  = dataFile.addVar ("pUvec_r",  ncDouble, vacuum_d);
      pwvecr_x.putVar (pUvec_r);
      
      NcVar psier_x = dataFile.addVar ("y_e_r",   ncDouble, ideal_d);
      psier_x.putVar (Psie_r);
      NcVar psiei_x = dataFile.addVar ("y_e_i",   ncDouble, ideal_d);
      psiei_x.putVar (Psie_i);
      NcVar zer_x   = dataFile.addVar ("Z_e_r",   ncDouble, ideal_d);
      zer_x.putVar (Ze_r);
      NcVar zei_x   = dataFile.addVar ("Z_e_i",   ncDouble, ideal_d);
      zei_x.putVar (Ze_i);

      NcVar ppsier_x = dataFile.addVar ("py_e_r",   ncDouble, ideal_d);
      ppsier_x.putVar (pPsie_r);
      NcVar ppsiei_x = dataFile.addVar ("py_e_i",   ncDouble, ideal_d);
      ppsiei_x.putVar (pPsie_i);
      NcVar pzer_x   = dataFile.addVar ("pZ_e_r",   ncDouble, ideal_d);
      pzer_x.putVar (pZe_r);
      NcVar pzei_x   = dataFile.addVar ("pZ_e_i",   ncDouble, ideal_d);
      pzei_x.putVar (pZe_i);
	  
      NcVar dW_x  = dataFile.addVar ("delta_W",   ncDouble, j_d);
      dW_x.putVar (deltaW);
      NcVar dWv_x = dataFile.addVar ("delta_W_p", ncDouble, j_d);
      dWv_x.putVar (deltaWp);
      NcVar dWp_x = dataFile.addVar ("delta_W_v", ncDouble, j_d);
      dWp_x.putVar (deltaWv);

      NcVar pdW_x  = dataFile.addVar ("pdelta_W",   ncDouble, j_d);
      pdW_x.putVar (pdeltaW);
      NcVar pdWv_x = dataFile.addVar ("pdelta_W_p", ncDouble, j_d);
      pdWv_x.putVar (pdeltaWp);
      NcVar pdWp_x = dataFile.addVar ("pdelta_W_v", ncDouble, j_d);
      pdWp_x.putVar (pdeltaWv);
      
      NcVar psiyr_x = dataFile.addVar ("y_surface_r", ncDouble, surface_d);
      psiyr_x.putVar (Psiy_r);
      NcVar psiyi_x = dataFile.addVar ("y_surface_i", ncDouble, surface_d);
      psiyi_x.putVar (Psiy_i);
      NcVar xyr_x   = dataFile.addVar ("Z_surface_r", ncDouble, surface_d);
      xyr_x.putVar (Xiy_r);
      NcVar xyi_x   = dataFile.addVar ("Z_surface_i", ncDouble, surface_d);
      xyi_x.putVar (Xiy_i);

      NcVar ppsiyr_x = dataFile.addVar ("py_surface_r", ncDouble, surface_d);
      ppsiyr_x.putVar (pPsiy_r);
      NcVar ppsiyi_x = dataFile.addVar ("py_surface_i", ncDouble, surface_d);
      ppsiyi_x.putVar (pPsiy_i);
      NcVar pxyr_x   = dataFile.addVar ("pZ_surface_r", ncDouble, surface_d);
      pxyr_x.putVar (pXiy_r);
      NcVar pxyi_x   = dataFile.addVar ("pZ_surface_i", ncDouble, surface_d);
      pxyi_x.putVar (pXiy_i);

      NcVar yar_x = dataFile.addVar ("ya_r", ncDouble, j_d);
      yar_x.putVar (ya_r);
      NcVar yai_x = dataFile.addVar ("ya_i", ncDouble, j_d);
      yai_x.putVar (ya_i);
      NcVar ybr_x = dataFile.addVar ("yb_r", ncDouble, j_d);
      ybr_x.putVar (yb_r);
      NcVar ybi_x = dataFile.addVar ("yb_i", ncDouble, j_d);
      ybi_x.putVar (yb_i);

      if (VIZ)
	{
	  NcVar ppvr_x = dataFile.addVar ("y_ideal_eig_r", ncDouble, v_d);
	  ppvr_x.putVar (PPV_r);
	  NcVar ppvi_x = dataFile.addVar ("y_ideal_eig_i", ncDouble, v_d);
	  ppvi_x.putVar (PPV_i);
	  NcVar zzvr_x = dataFile.addVar ("Z_ideal_eig_r", ncDouble, v_d);
	  zzvr_x.putVar (ZZV_r);
	  NcVar zzvi_x = dataFile.addVar ("Z_ideal_eig_i", ncDouble, v_d);
	  zzvi_x.putVar (ZZV_i);

	  NcVar pppvr_x = dataFile.addVar ("py_ideal_eig_r", ncDouble, v_d);
	  pppvr_x.putVar (pPPV_r);
	  NcVar pppvi_x = dataFile.addVar ("py_ideal_eig_i", ncDouble, v_d);
	  pppvi_x.putVar (pPPV_i);
	  NcVar pzzvr_x = dataFile.addVar ("pZ_ideal_eig_r", ncDouble, v_d);
	  pzzvr_x.putVar (pZZV_r);
	  NcVar pzzvi_x = dataFile.addVar ("pZ_ideal_eig_i", ncDouble, v_d);
	  pzzvi_x.putVar (pZZV_i);
	}
   
      dataFile.close ();
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Vertical/Virtical.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

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
  delete[] iGmat_r; delete[] iGmat_i; delete[] Iant_r;  delete[] Iant_i;

  delete[] Ammp_r; delete[] Bmmp_r; delete[] Cmmp_r; delete[] Dmmp_r;
  delete[] Ammp_i; delete[] Bmmp_i; delete[] Cmmp_i; delete[] Dmmp_i;
  delete[] Atest;  delete[] BCtest; delete[] Dtest;

  delete[] PPPsi_r; delete[] PPPsi_i; delete[] ZZZ_r; delete[] ZZZ_i;

  delete[] Umat_r; delete[] Umat_i; delete[] Uant_r; delete[] Uant_i;
  delete[] Psie_r; delete[] Psie_i; delete[] Ze_r;   delete[] Ze_i;
  delete[] Psiy_r; delete[] Psiy_i; delete[] Xiy_r;  delete[] Xiy_i;

  delete[] pUmat_r; delete[] pUmat_i; delete[] pUant_r; delete[] pUant_i;
  delete[] pPsie_r; delete[] pPsie_i; delete[] pZe_r;   delete[] pZe_i;
  delete[] pPsiy_r; delete[] pPsiy_i; delete[] pXiy_r;  delete[] pXiy_i;

  delete[] PPV_r; delete[] PPV_i; delete[] ZZV_r; delete ZZV_i;

  delete[] pPPV_r; delete[] pPPV_i; delete[] pZZV_r; delete pZZV_i;

  delete[] ya_r; delete[] ya_i; delete[] yb_r; delete[] yb_i;
}
