// Stage2.cpp

// PROGRAM ORGANIATION:
//
// void Flux:: Stage2             ()
// void Flux:: WriteStage2Netcdfc ()

#include "Flux.h"

// ####################################################
// Function to input Stage1 data and output Stage2 data
// ####################################################
void Flux::Stage2 ()
{
  // ................................
  // Read data for Stage2 calculation
  // ................................
  Stage2ReadData ();
  fflush (stdout);

  // ..........................
  // Calculate Stage2 q profile
  // ..........................
  Stage2CalcQ ();
  fflush (stdout);

  // ....................................
  // Calculate Stage2 straight angle data
  // ....................................
  Stage2CalcStraightAngle ();
  fflush (stdout);

  // .......................
  // Calculate boundary data
  // .......................
  Stage2CalcBoundary ();
  fflush (stdout);
  
  // .........................
  // Output Stage2 NETCDF file
  // .........................
  WriteStage2Netcdfc ();

  // ...........
  // Free memory
  // ...........
  delete[] RPTS;  delete[] ZPTS;
  delete[] RBPTS; delete[] ZBPTS;
  delete[] RLPTS; delete[] ZLPTS;

  delete[] PSIN; delete[] G;   delete[] Pr;
  delete[] GGp;  delete[] Prp; delete[] Q;

  
  delete[] s; delete[] Rs; delete[] PSIAXIS;

  delete[] P;    delete[] RP;  delete[] rP;   delete[] GP;
  delete[] QGP;  delete[] QP;  delete[] PP;   delete[] GPP;
  delete[] PPP;  delete[] S;   delete[] QX;   delete[] PPX; 
  delete[] PsiN; delete[] FP;  delete[] GPX;  delete[] th; 

  gsl_matrix_free (RRst); gsl_matrix_free (ZZst);
  gsl_matrix_free (RRr);  gsl_matrix_free (RRth);
  gsl_matrix_free (ZZr);  gsl_matrix_free (ZZth);
  gsl_matrix_free (Jac);  gsl_matrix_free (Jax);

  delete[] Rbndry; delete[] Zbndry; delete[] dRdtheta; delete[] dZdtheta; 
}

// ####################################
// Function to write Stage2 NETCDF file
// ####################################
void Flux::WriteStage2Netcdfc ()
{
  try
    {
      NcFile dataFile ("Outputs/Flux/Stage2.nc", NcFile::replace);

      double parameters[10];
      parameters[0] = R0;
      parameters[1] = B0;
      parameters[2] = RLEFT;
      parameters[3] = ZLOW;
      parameters[4] = RRIGHT;
      parameters[5] = ZHIGH;
      parameters[6] = Raxis;
      parameters[7] = Zaxis;
      parameters[8] = Psic;
      parameters[9] = Rbound;

      int ipara[4];
      ipara[0] = ia;
      ipara[1] = ic;
      ipara[2] = jc;
      ipara[3] = L;

      NcDim para_d = dataFile.addDim ("r_index", 10);
      NcVar p_x    = dataFile.addVar ("para",  ncDouble, para_d);
      p_x.putVar (parameters);

      NcDim ipara_d = dataFile.addDim ("i_index", 4);
      NcVar ip_x    = dataFile.addVar ("i_para", ncInt, ipara_d);
      ip_x.putVar (ipara);

      NcDim bound_d   = dataFile.addDim ("i_bound", NBPTS);
      NcVar bound_r_x = dataFile.addVar ("RBPTS", ncDouble, bound_d);
      bound_r_x.putVar (RBPTS);
      NcVar bound_z_x = dataFile.addVar ("ZBPTS", ncDouble, bound_d);
      bound_z_x.putVar (ZBPTS);

      NcDim R_d = dataFile.addDim ("i", NRPTS);
      NcDim Z_d = dataFile.addDim ("j", NZPTS);
      NcVar r_x = dataFile.addVar ("R", ncDouble, R_d);
      r_x.putVar (RPTS);
      NcVar z_x = dataFile.addVar ("Z", ncDouble, Z_d);
      z_x.putVar (ZPTS);
      
      int cnt = 0;
      double* DATA = new double[NRPTS*NZPTS];
      for (int i = 0; i < NRPTS; i++)
	for (int j = 0; j < NZPTS; j++)
	  {
	    DATA[cnt] = PSIARRAY (i, j);
	    cnt++;
	  }

      vector<NcDim> psi_d;
      psi_d.push_back (R_d);
      psi_d.push_back (Z_d);
      NcVar psi_x = dataFile.addVar ("PSI", ncDouble, psi_d);
      psi_x.putVar (DATA);

      cnt = 0;
      double* DATA1 = new double[NRPTS*NZPTS];
      for (int i = 0; i < NRPTS; i++)
	for (int j = 0; j < NZPTS; j++)
	  {
	    DATA1[cnt] = GetPsiR (RPTS[i], ZPTS[j]);
	    cnt++;
	  }

      NcVar psir_x = dataFile.addVar ("PSI_R", ncDouble, psi_d);
      psir_x.putVar (DATA1);

      cnt = 0;
      double* DATA2 = new double[NRPTS*NZPTS];
      for (int i = 0; i < NRPTS; i++)
	for (int j = 0; j < NZPTS; j++)
	  {
	    DATA2[cnt] = GetPsiZ (RPTS[i], ZPTS[j]);
	    cnt++;
	  }

      NcVar psiz_x = dataFile.addVar ("PSI_Z", ncDouble, psi_d);
      psiz_x.putVar (DATA2);

      NcDim L_d = dataFile.addDim ("l", L);
      NcVar s_x = dataFile.addVar ("s",   ncDouble, L_d);
      s_x.putVar (s);
      NcVar Rs_x = dataFile.addVar ("Rs", ncDouble, L_d);
      Rs_x.putVar (Rs);

      NcDim Q_d = dataFile.addDim ("N_PSI", NPSI);
      NcVar PsiN_x = dataFile.addVar  ("PsiN",    ncDouble, Q_d);
      PsiN_x.putVar (PsiN);

      NcVar q_x = dataFile.addVar     ("q",       ncDouble, Q_d);
      q_x.putVar (QP);

      NcVar qx_x = dataFile.addVar    ("q_x",     ncDouble, Q_d);
      qx_x.putVar (QX);

      NcVar g_x = dataFile.addVar     ("g",       ncDouble, Q_d);
      g_x.putVar (GP);
      
      NcVar f_x = dataFile.addVar     ("f",       ncDouble, Q_d);
      f_x.putVar (FP);

      NcVar pp_x = dataFile.addVar    ("P",       ncDouble, Q_d);
      pp_x.putVar (PP);

      NcVar gpsi_x = dataFile.addVar  ("g_psi",   ncDouble, Q_d);
      gpsi_x.putVar (GPP);

      NcVar ppsi_x = dataFile.addVar  ("P_psi",   ncDouble, Q_d);
      ppsi_x.putVar (PPP);

      NcVar gpsix_x = dataFile.addVar ("g_psi_x", ncDouble, Q_d);
      gpsi_x.putVar (GPX);

      NcVar ppsix_x = dataFile.addVar ("P_psi_x", ncDouble, Q_d);
      ppsi_x.putVar (PPX);

      double* rra = new double[NPSI];
      for (int j = 0; j < NPSI; j++)
	rra [j] = rP[j] /ra;
      NcVar rr_x = dataFile.addVar ("r", ncDouble, Q_d);
      rr_x.putVar (rra);

      NcDim T_d = dataFile.addDim ("k", NTHETA);
      double* tth = new double[NTHETA];
      for (int j = 0; j < NTHETA; j++)
	tth[j] = th[j] /M_PI;
      NcVar th_x = dataFile.addVar ("theta", ncDouble, T_d);
      th_x.putVar (tth);

      cnt = 0;
      double* DATA5 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA5[cnt] = gsl_matrix_get (RRst, i, j);
	    cnt++;
	  }

      vector<NcDim> rrst_d;
      rrst_d.push_back (Q_d);
      rrst_d.push_back (T_d);
      NcVar rrst_x = dataFile.addVar ("RRst", ncDouble, rrst_d);
      rrst_x.putVar (DATA5);

      cnt = 0;
      double* DATA6 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA6[cnt] = gsl_matrix_get (ZZst, i, j);
	    cnt++;
	  }
      NcVar zzst_x = dataFile.addVar ("ZZst", ncDouble, rrst_d);
      zzst_x.putVar (DATA6);

      cnt = 0;
      double* DATA7 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA7[cnt] = gsl_matrix_get (RRr, i, j);
	    cnt++;
	  }
      NcVar drdr_x = dataFile.addVar ("dRdr", ncDouble, rrst_d);
      drdr_x.putVar (DATA7);

      cnt = 0;
      double* DATA8 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA8[cnt] = gsl_matrix_get (RRth, i, j);
	    cnt++;
	  }
      NcVar drdth_x = dataFile.addVar ("dRdth", ncDouble, rrst_d);
      drdth_x.putVar (DATA8);

      cnt = 0;
      double* DATA9 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA9[cnt] = gsl_matrix_get (ZZr, i, j);
	    cnt++;
	  }
      NcVar dzdr_x = dataFile.addVar ("dZdr", ncDouble, rrst_d);
      dzdr_x.putVar (DATA9);

      cnt = 0;
      double* DATA10 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA10[cnt] = gsl_matrix_get (ZZth, i, j);
	    cnt++; 
	  }
      NcVar dzdth_x = dataFile.addVar ("dZdth", ncDouble, rrst_d);
      dzdth_x.putVar (DATA10);

      cnt = 0;
      double* DATA11 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA11[cnt] = gsl_matrix_get (Jac, i, j);
	    cnt++;
	  }
      NcVar jac_x = dataFile.addVar ("Jac", ncDouble, rrst_d);
      jac_x.putVar (DATA11);

      cnt = 0;
      double* DATA12 = new double[NPSI*NTHETA];
      for (int i = 0; i < NPSI; i++)
	for (int j = 0; j < NTHETA; j++)
	  {
	    DATA12[cnt] = gsl_matrix_get (Jax, i, j);
	    cnt++;
	  }
       NcVar jax_x = dataFile.addVar ("Jax", ncDouble, rrst_d);
       jax_x.putVar (DATA12);

       NcVar tbound_x   = dataFile.addVar ("tbound",   ncDouble, T_d);
       tbound_x.putVar (th);
       
       NcVar rbound_x   = dataFile.addVar ("Rbound",   ncDouble, T_d);
       rbound_x.putVar (Rbndry);

       NcVar zbound_x   = dataFile.addVar ("Zbound",   ncDouble, T_d);
       zbound_x.putVar (Zbndry);

       NcVar drdtheta_x = dataFile.addVar ("dRdtheta", ncDouble, T_d);
       drdtheta_x.putVar (dRdtheta);

       NcVar dzdtheta_x = dataFile.addVar ("dZdtheta", ncDouble, T_d);
       dzdtheta_x.putVar (dZdtheta);

       delete[] rra;    delete[] tth;     
       delete[] DATA;   delete[] DATA1,  delete[] DATA2;   
       delete[] DATA5;  delete[] DATA6;  delete[] DATA7;
       delete[] DATA8;  delete[] DATA9;  delete[] DATA10;
       delete[] DATA11; delete[] DATA12;
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Flux/Stage2.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }
}

