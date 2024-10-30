// Netcdf.cpp

#include "Equilibrium.h"

// #################################################
// Function to write equilibrium data to netcdf file
// #################################################
void Equilibrium::WriteNetcdf (double sa)
{
  printf ("Writing equilibrium data to netcdf file Outputs/Equilibrium/Equilibrium.nc:\n");

  double para[2], Input[15], Beta[4];

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
  
  para[0] = epsa;
  para[1] = sa;

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
  
      NcDim i_d = dataFile.addDim ("Ni", 15);
      NcDim p_d = dataFile.addDim ("Np", 2);
      NcDim b_d = dataFile.addDim ("Nb", 4);
      NcDim r_d = dataFile.addDim ("Nr", Nr+1);
      NcDim s_d = dataFile.addDim ("Ns", Ns+1);
      NcDim f_d = dataFile.addDim ("Nf", Nf);
      NcDim w_d = dataFile.addDim ("Nw", Nw+1);
 
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
 
      NcVar Hn_x  = dataFile.addVar ("Hn",  ncDouble, shape_d);
      Hn_x.putVar (HHfunc.data());
      NcVar Hnp_x = dataFile.addVar ("Hnp", ncDouble, shape_d);
      Hnp_x.putVar (HPfunc.data());
      NcVar Vn_x  = dataFile.addVar ("Vn",  ncDouble, shape_d);
      Vn_x.putVar (VVfunc.data());
      NcVar Vnp_x = dataFile.addVar ("Vnp", ncDouble, shape_d);
      Vnp_x.putVar (VPfunc.data());

      NcVar R_x  = dataFile.addVar ("R",     ncDouble, flux_d);
      R_x.putVar (RR.data());
      NcVar Z_x  = dataFile.addVar ("Z",     ncDouble, flux_d);
      Z_x.putVar (ZZ.data());
      NcVar Rw_x = dataFile.addVar ("Rw",    ncDouble, flux_d);
      Rw_x.putVar (RRw.data());
      NcVar Zw_x = dataFile.addVar ("Zw",    ncDouble, flux_d);
      Zw_x.putVar (ZZw.data());
      NcVar rr_x = dataFile.addVar ("rr",    ncDouble, flux_d);
      rr_x.putVar (rvals.data());
      NcVar t_x  = dataFile.addVar ("theta", ncDouble, flux_d);
      t_x.putVar (thvals.data());
      NcVar w_x  = dataFile.addVar ("omega", ncDouble, flux_d);
      w_x.putVar (wvals.data());

      NcVar n_x   = dataFile.addVar ("n",   ncDouble, s_d);
      n_x.putVar (npol);
      NcVar Hna_x = dataFile.addVar ("Hna", ncDouble, s_d);
      Hna_x.putVar (Hna);
      NcVar Vna_x = dataFile.addVar ("Vna", ncDouble, s_d);
      Vna_x.putVar (Vna);

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
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Equilibrium/Equilibrium.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] npol; delete[] Hna; delete[] Vna;
}
