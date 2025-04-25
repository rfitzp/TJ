// Netcdf.cpp

#include "Equilibrium.h"

#define NINPUT 18
#define NPARA 3
#define NBETA 4

// #################################################
// Function to write equilibrium data to netcdf file
// #################################################
void Equilibrium::WriteNetcdf (double sa)
{
  printf ("Writing equilibrium data to netcdf file Outputs/Equilibrium/Equilibrium.nc:\n");

  double para[NPARA], Input[NINPUT], Beta[NBETA];

  double* Hna   = new double[Ns+1];
  double* Vna   = new double[Ns+1];
  double* npol  = new double[Ns+1];

  Input[0]  = qc;
  Input[1]  = qa;
  Input[2]  = nu;
  Input[3]  = pc;
  Input[4]  = mu;
  Input[5]  = epsa;
  Input[6]  = eps;
  Input[7]  = double (Ns);
  Input[8]  = double (Nr);
  Input[9]  = double (Nf);
  Input[10] = double (Nw);
  Input[11] = acc;
  Input[12] = h0;
  Input[13] = hmin;
  Input[14] = hmax;
  Input[15] = double (SRC);
  Input[16] = bw;
  Input[17] = tilt;
  
  para[0] = epsa;
  para[1] = sa;
  para[2] = bw;

  Beta[0] = li;
  Beta[1] = betat;
  Beta[2] = betap;
  Beta[3] = betaN;
    
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
    npol[n] = double (n);
 
  try
    {
      NcFile dataFile ("../Outputs/Equilibrium/Equilibrium.nc", NcFile::replace);

      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);
  
      NcDim i_d = dataFile.addDim ("Ni",    NINPUT);
      NcDim p_d = dataFile.addDim ("Np",    NPARA);
      NcDim b_d = dataFile.addDim ("Nb",    NBETA);
      NcDim r_d = dataFile.addDim ("Nr",    Nr+1);
      NcDim s_d = dataFile.addDim ("Ns",    Ns+1);
      NcDim f_d = dataFile.addDim ("Nf",    Nf);
      NcDim w_d = dataFile.addDim ("Nw",    Nw+1);
      NcDim c_d = dataFile.addDim ("ncoil", ncoil);
      NcDim e_d = dataFile.addDim ("Neq",   2*Nf);
 
      vector<NcDim> shape_d;
      shape_d.push_back (s_d);
      shape_d.push_back (r_d);

      vector<NcDim> flux_d;
      flux_d.push_back (f_d);
      flux_d.push_back (w_d);

      NcVar i_x = dataFile.addVar ("InputParameters", ncDouble, i_d);
      i_x.putVar (Input);
      NcVar p_x = dataFile.addVar ("para",            ncDouble, p_d);
      p_x.putVar (para);
      NcVar b_x = dataFile.addVar ("beta",            ncDouble, b_d);
      b_x.putVar (Beta);

      NcVar r_x   = dataFile.addVar ("r",    ncDouble, r_d);
      r_x.putVar (rr);
      NcVar psi_x = dataFile.addVar ("Psi",  ncDouble, r_d);
      psi_x.putVar (Psi);
      NcVar psn_x = dataFile.addVar ("PsiN", ncDouble, r_d);
      psn_x.putVar (PsiN);
      NcVar g2_x  = dataFile.addVar ("g_2",  ncDouble, r_d);
      g2_x.putVar (g2);
      NcVar p2_x  = dataFile.addVar ("p_2",  ncDouble, r_d);
      p2_x.putVar (p2);
      NcVar pp_x  = dataFile.addVar ("pp",   ncDouble, r_d);
      pp_x.putVar (pp);
      NcVar ppp_x = dataFile.addVar ("ppp",  ncDouble, r_d);
      ppp_x.putVar (ppp);
      NcVar f1_x  = dataFile.addVar ("f_1",  ncDouble, r_d);
      f1_x.putVar (f1);
      NcVar f3_x  = dataFile.addVar ("f_3",  ncDouble, r_d);
      f3_x.putVar (f3);
      NcVar f_x   = dataFile.addVar ("f",    ncDouble, r_d);
      f_x.putVar (ff);
      NcVar q0_x  = dataFile.addVar ("q_0",  ncDouble, r_d);
      q0_x.putVar (q0);
      NcVar q2_x  = dataFile.addVar ("q_2",  ncDouble, r_d);
      q2_x.putVar (q2);
      NcVar It_x  = dataFile.addVar ("I_t",  ncDouble, r_d);
      It_x.putVar (It);
      NcVar Ip_x  = dataFile.addVar ("I_p",  ncDouble, r_d);
      Ip_x.putVar (Ip);
      NcVar Jt_x  = dataFile.addVar ("J_t",  ncDouble, r_d);
      Jt_x.putVar (Jt);
      NcVar Jp_x  = dataFile.addVar ("J_p",  ncDouble, r_d);
      Jp_x.putVar (Jp);
      NcVar q_x   = dataFile.addVar ("q",    ncDouble, r_d);
      q_x.putVar (q2);
      NcVar qq_x  = dataFile.addVar ("qq",   ncDouble, r_d);
      qq_x.putVar (qq);
      NcVar qqq_x = dataFile.addVar ("qqq",  ncDouble, r_d);
      qqq_x.putVar (qqq);
      NcVar s_x   = dataFile.addVar ("s",    ncDouble, r_d);
      s_x.putVar (s);
      NcVar s2_x  = dataFile.addVar ("s2",   ncDouble, r_d);
      s2_x.putVar (s2);
      NcVar s0_x  = dataFile.addVar ("s0",   ncDouble, r_d);
      s0_x.putVar (s0);
      NcVar S1_x  = dataFile.addVar ("S1",   ncDouble, r_d);
      S1_x.putVar (S1);
      NcVar S2_x  = dataFile.addVar ("S2",   ncDouble, r_d);
      S2_x.putVar (S2);
      NcVar S3_x  = dataFile.addVar ("S3",   ncDouble, r_d);
      S3_x.putVar (S3);
      NcVar P1_x  = dataFile.addVar ("P1",   ncDouble, r_d);
      P1_x.putVar (P1);
      NcVar P1a_x = dataFile.addVar ("P1a",  ncDouble, r_d);
      P1a_x.putVar (P1a);
      NcVar P2_x  = dataFile.addVar ("P2",   ncDouble, r_d);
      P2_x.putVar (P2);
      NcVar P2a_x = dataFile.addVar ("P2a",  ncDouble, r_d);
      P2a_x.putVar (P2a);
      NcVar P3_x  = dataFile.addVar ("P3",   ncDouble, r_d);
      P3_x.putVar (P3);
      NcVar P3a_x = dataFile.addVar ("P3a",  ncDouble, r_d);
      P3a_x.putVar (P3a);
      NcVar T_x   = dataFile.addVar ("T",    ncDouble, r_d);
      T_x.putVar (Tf);
      NcVar mP_x  = dataFile.addVar ("mu0P", ncDouble, r_d);
      mP_x.putVar (mu0P);
      NcVar DI_x  = dataFile.addVar ("DI",   ncDouble, r_d);
      DI_x.putVar (DI);
      NcVar DR_x  = dataFile.addVar ("DR",   ncDouble, r_d);
      DR_x.putVar (DR);
      NcVar ne_x  = dataFile.addVar ("ne",   ncDouble, r_d);
      ne_x.putVar (ne);
      NcVar Te_x  = dataFile.addVar ("Te",   ncDouble, r_d);
      Te_x.putVar (Te);
      NcVar nep_x = dataFile.addVar ("nep",  ncDouble, r_d);
      nep_x.putVar (nep);
      NcVar Tep_x = dataFile.addVar ("Tep",  ncDouble, r_d);
      Tep_x.putVar (Tep);
 
      NcVar Hn_x  = dataFile.addVar ("Hn",  ncDouble, shape_d);
      Hn_x.putVar (HHfunc.data());
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (HPfunc.data());
      NcVar Vn_x  = dataFile.addVar ("Vn",  ncDouble, shape_d);
      Vn_x.putVar (VVfunc.data());
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (VPfunc.data());

      NcVar n_x   = dataFile.addVar ("n",   ncDouble, s_d);
      n_x.putVar (npol);
      NcVar Hna_x = dataFile.addVar ("Hna", ncDouble, s_d);
      Hna_x.putVar (Hna);
      NcVar Vna_x = dataFile.addVar ("Vna", ncDouble, s_d);
      Vna_x.putVar (Vna);

      if (VIZ)
	{
	  NcVar R_x    = dataFile.addVar ("R",     ncDouble, flux_d);
	  R_x.putVar (RR.data());
	  NcVar Z_x    = dataFile.addVar ("Z",     ncDouble, flux_d);
	  Z_x.putVar (ZZ.data());
	  NcVar dRdr_x = dataFile.addVar ("dRdr",  ncDouble, flux_d);
	  dRdr_x.putVar (dRdr.data());
	  NcVar dRdt_x = dataFile.addVar ("dRdt",  ncDouble, flux_d);
	  dRdt_x.putVar (dRdt.data());
	  NcVar dZdr_x = dataFile.addVar ("dZdr",  ncDouble, flux_d);
	  dZdr_x.putVar (dZdr.data());
	  NcVar dZdt_x = dataFile.addVar ("dZdt",  ncDouble, flux_d);
	  dZdt_x.putVar (dZdt.data());
	  NcVar Jac_x  = dataFile.addVar ("Jac",   ncDouble, flux_d);
	  Jac_x.putVar (Jac.data());
	  NcVar Jax_x  = dataFile.addVar ("Jax",   ncDouble, flux_d);
	  Jax_x.putVar (Jax.data());
	  NcVar Rw_x   = dataFile.addVar ("Rw",    ncDouble, flux_d);
	  Rw_x.putVar (RRw.data());
	  NcVar Zw_x   = dataFile.addVar ("Zw",    ncDouble, flux_d);
	  Zw_x.putVar (ZZw.data());
	  NcVar rr_x   = dataFile.addVar ("rr",    ncDouble, flux_d);
	  rr_x.putVar (rvals.data());
	  NcVar t_x    = dataFile.addVar ("theta", ncDouble, flux_d);
	  t_x.putVar (thvals.data());
	  NcVar w_x    = dataFile.addVar ("omega", ncDouble, flux_d);
	  w_x.putVar (wvals.data());

	  NcVar req_x  = dataFile.addVar ("r_eq",     ncDouble, e_d);
	  req_x.putVar (req);
	  NcVar weq_x  = dataFile.addVar ("omega_eq", ncDouble, e_d);
	  weq_x.putVar (weq);
	  NcVar teq_x  = dataFile.addVar ("theta_eq", ncDouble, e_d);
	  teq_x.putVar (teq);
	  NcVar Req_x  = dataFile.addVar ("R_eq",     ncDouble, e_d);
	  Req_x.putVar (Req);
	  NcVar Zeq_x  = dataFile.addVar ("Z_eq",     ncDouble, e_d);
	  Zeq_x.putVar (Zeq);
	  NcVar BReq_x = dataFile.addVar ("BR_eq",    ncDouble, e_d);
	  BReq_x.putVar (BReq);
	  NcVar neeq_x = dataFile.addVar ("ne_eq",    ncDouble, e_d);
	  neeq_x.putVar (neeq);
	  NcVar Teeq_x = dataFile.addVar ("Te_eq",    ncDouble, e_d);
	  Teeq_x.putVar (Teeq);
	  NcVar Rreq_x = dataFile.addVar ("dRdr_eq",  ncDouble, e_d);
	  Rreq_x.putVar (dRdreq);
	  NcVar Rteq_x = dataFile.addVar ("dRdt_eq",  ncDouble, e_d);
	  Rteq_x.putVar (dRdteq);
	  NcVar Zreq_x = dataFile.addVar ("dZdr_eq",  ncDouble, e_d);
	  Zreq_x.putVar (dZdreq);
	  NcVar Zteq_x = dataFile.addVar ("dZdt_eq",  ncDouble, e_d);
	  Zteq_x.putVar (dZdteq);
	}

      NcVar Rbound_x   = dataFile.addVar ("Rbound",    ncDouble, w_d);
      Rbound_x.putVar (Rbound);
      NcVar Zbound_x   = dataFile.addVar ("Zbound",    ncDouble, w_d);
      Zbound_x.putVar (Zbound);
      NcVar tbound_x   = dataFile.addVar ("tbound",    ncDouble, w_d);
      tbound_x.putVar (tbound);
      NcVar wbound_x   = dataFile.addVar ("wbound",    ncDouble, w_d);
      wbound_x.putVar (wbound);
      NcVar R2b_x      = dataFile.addVar ("R2bound",   ncDouble, w_d);
      R2b_x.putVar (R2b);
      NcVar grr2b_x    = dataFile.addVar ("grr2bound", ncDouble, w_d);
      grr2b_x.putVar (grr2b);
      NcVar dRdtheta_x = dataFile.addVar ("dRdtheta",  ncDouble, w_d);
      dRdtheta_x.putVar (dRdtheta);
      NcVar dZdtheta_x = dataFile.addVar ("dZdtheta",  ncDouble, w_d);
      dZdtheta_x.putVar (dZdtheta);

      NcVar Rwall_x = dataFile.addVar ("Rwall", ncDouble, w_d);
      Rwall_x.putVar (Rwall);
      NcVar Zwall_x = dataFile.addVar ("Zwall", ncDouble, w_d);
      Zwall_x.putVar (Zwall);
      NcVar wwall_x = dataFile.addVar ("wwall", ncDouble, w_d);
      wwall_x.putVar (wwall);
 
      NcVar Rcoil_x = dataFile.addVar ("Rcoil", ncDouble, c_d);
      Rcoil_x.putVar (Rcoil);
      NcVar Zcoil_x = dataFile.addVar ("Zcoil", ncDouble, c_d);
      Zcoil_x.putVar (Zcoil);
      NcVar Icoil_x = dataFile.addVar ("Icoil", ncDouble, c_d);
      Icoil_x.putVar (Icoil);

      dataFile.close ();
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Equilibrium/Equilibrium.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] npol; delete[] Hna; delete[] Vna;
}
