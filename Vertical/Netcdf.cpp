// Netcdf.cpp

#include "Vertical.h"

#define NINPUT 21

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
      NcVar S5_x  = dataFile.getVar ("S5");
      NcVar P4_x  = dataFile.getVar ("P3a");
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
      S5   = new double[Nr+1];
      Sig  = new double[Nr+1];
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
      S5_x .getVar (S5);
      P4_x .getVar (Sig);
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

  double Input[NINPUT];

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
  
  try
    {
      NcFile dataFile ("../Outputs/Vertical/Vertical.nc", NcFile::replace);

      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim i_d = dataFile.addDim ("Ni", NINPUT);
      NcDim r_d = dataFile.addDim ("Nr", Nr+1);
      NcDim s_d = dataFile.addDim ("Ns", Ns+1);
      NcDim j_d = dataFile.addDim ("J",  J);

      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

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

      dataFile.close ();
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Vertical/Virtical.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }
}
