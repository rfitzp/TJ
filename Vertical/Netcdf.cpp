// Netcdf.cpp

#include "Vertical.h"

#define NINPUT 27

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
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P3_x  = dataFile.getVar ("P3");
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
      P1   = new double[Nr+1];
      P2   = new double[Nr+1];
      P3   = new double[Nr+1];
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
      P1_x .getVar (P1);
      P2_x .getVar (P2);
      P3_x .getVar (P3);
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

	  NcVar dRdr_x = dataFile.getVar ("dRdr");
	  NcVar dRdt_x = dataFile.getVar ("dRdt");
	  NcVar dZdr_x = dataFile.getVar ("dZdr");
	  NcVar dZdt_x = dataFile.getVar ("dZdt");
	  
	  Nf = f_d.getSize ();
	  
	  RR.resize     (Nf, Nw+1);
	  ZZ.resize     (Nf, Nw+1);
	  rvals.resize  (Nf, Nw+1);
	  thvals.resize (Nf, Nw+1);

	  dRdr.resize (Nf, Nw+1);
	  dRdt.resize (Nf, Nw+1);
	  dZdr.resize (Nf, Nw+1);
	  dZdt.resize (Nf, Nw+1);
	  
	  double* RRdata = new double[Nf*(Nw+1)];
	  double* ZZdata = new double[Nf*(Nw+1)];
	  double* rrdata = new double[Nf*(Nw+1)];
	  double* ttdata = new double[Nf*(Nw+1)];

	  double* dRdrdata = new double[Nf*(Nw+1)];
	  double* dRdtdata = new double[Nf*(Nw+1)];
	  double* dZdrdata = new double[Nf*(Nw+1)];
	  double* dZdtdata = new double[Nf*(Nw+1)];	  
	  
	  RR_x.getVar (RRdata);
	  ZZ_x.getVar (ZZdata);
	  rr_x.getVar (rrdata);
	  tt_x.getVar (ttdata);

	  dRdr_x.getVar (dRdrdata);
	  dRdt_x.getVar (dRdtdata);
	  dZdr_x.getVar (dZdrdata);
	  dZdt_x.getVar (dZdtdata);
	  
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

		dRdr(n, i) = dRdrdata[i + n*(Nw+1)];
		dRdt(n, i) = dRdtdata[i + n*(Nw+1)];
		dZdr(n, i) = dZdrdata[i + n*(Nw+1)];
		dZdt(n, i) = dZdtdata[i + n*(Nw+1)];
	      }

	  delete[] RRdata;   delete[] ZZdata;   delete[] rrdata;   delete[] ttdata;
	  delete[] dRdrdata; delete[] dRdtdata; delete[] dZdrdata; delete[] dZdtdata;

	  NcVar req_x    = dataFile.getVar ("r_eq");
	  NcVar teq_x    = dataFile.getVar ("theta_eq");
	  NcVar Req_x    = dataFile.getVar ("R_eq");
	  NcVar Zeq_x    = dataFile.getVar ("Z_eq");
	  NcVar BReq_x   = dataFile.getVar ("BR_eq");
	  NcVar neeq_x   = dataFile.getVar ("ne_eq");
	  NcVar Teeq_x   = dataFile.getVar ("Te_eq");
	  NcVar dRdreq_x = dataFile.getVar ("dRdr_eq");
	  NcVar dRdteq_x = dataFile.getVar ("dRdt_eq");
	  NcVar dZdreq_x = dataFile.getVar ("dZdr_eq");
	  NcVar dZdteq_x = dataFile.getVar ("dZdt_eq");

	  req    = new double[2*Nf];
	  teq    = new double[2*Nf];
	  Req    = new double[2*Nf];
	  Zeq    = new double[2*Nf];
	  BReq   = new double[2*Nf];
	  neeq   = new double[2*Nf];
	  Teeq   = new double[2*Nf];
	  dRdreq = new double[2*Nf];
	  dRdteq = new double[2*Nf];
	  dZdreq = new double[2*Nf];
	  dZdteq = new double[2*Nf];

	  req_x   .getVar (req);
	  teq_x   .getVar (teq);
	  Req_x   .getVar (Req);
	  Zeq_x   .getVar (Zeq);
	  BReq_x  .getVar (BReq);
	  neeq_x  .getVar (neeq);
	  Teeq_x  .getVar (Teeq);
	  dRdreq_x.getVar (dRdreq);
	  dRdteq_x.getVar (dRdteq);
	  dZdreq_x.getVar (dZdreq);
	  dZdteq_x.getVar (dZdteq);

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

