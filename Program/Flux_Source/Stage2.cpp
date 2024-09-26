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
  //Stage2CalcQ ();
  fflush (stdout);

  // ....................................
  // Calculate Stage2 straight angle data
  // ....................................
  // Stage2CalcStraightAngle ();
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
}

// ####################################
// Function to write Stage2 NETCDF file
// ####################################
void Flux::WriteStage2Netcdfc ()
{
  try
    {
      NcFile dataFile ("Outputs/Stage2/Stage2.nc", NcFile::replace);

      double parameters[8];
      parameters[0] = R0;
      parameters[1] = B0;
      parameters[2] = RLEFT;
      parameters[3] = ZLOW;
      parameters[4] = RRIGHT;
      parameters[5] = ZHIGH;
      parameters[6] = Raxis;
      parameters[7] = Zaxis;

      NcDim para_d = dataFile.addDim ("index", 8);
      NcVar p_x    = dataFile.addVar ("para",  ncDouble, para_d);
      p_x.putVar (parameters);

      NcDim bound_d   = dataFile.addDim ("i_bound", NBPTS);
      NcVar bound_r_x = dataFile.addVar ("RBPTS", ncDouble, bound_d);
      bound_r_x.putVar (RBPTS);
      NcVar bound_z_x = dataFile.addVar ("ZBPTS", ncDouble, bound_d);
      bound_z_x.putVar (ZBPTS);

      NcDim R_d   = dataFile.addDim ("i", NRPTS);
      NcDim Z_d   = dataFile.addDim ("j", NRPTS);
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
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Plots/Equilbrium.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }
 
  /*
  
 
  // Psi_N
  int Q_d, PN;
  err += nc_def_dim (dataFile, "N_PSI", NPSI, &Q_d);
  err += nc_def_var (dataFile, "PsiN", NC_DOUBLE, 1, &Q_d, &PN);

  // q
  int Q;
  err += nc_def_var (dataFile, "q", NC_DOUBLE, 1, &Q_d, &Q);

  // gamma
  int GMp;
  err += nc_def_var (dataFile, "gamma", NC_DOUBLE, 1, &Q_d, &GMp);

  // q_x
  int Qx;
  err += nc_def_var (dataFile, "q_x", NC_DOUBLE, 1, &Q_d, &Qx);
 
  // g
  int G;
  err += nc_def_var (dataFile, "g", NC_DOUBLE, 1, &Q_d, &G);

  // f
  int F;
  err += nc_def_var (dataFile, "f", NC_DOUBLE, 1, &Q_d, &F);

  // P
  int Px;
  err += nc_def_var (dataFile, "P", NC_DOUBLE, 1, &Q_d, &Px);
 
  // g_psi
  int Gp;
  err += nc_def_var (dataFile, "g_psi", NC_DOUBLE, 1, &Q_d, &Gp);
  
  // P_psi
  int Pp;
  err += nc_def_var (dataFile, "P_psi", NC_DOUBLE, 1, &Q_d, &Pp);

  // g_psi_x
  int Gpx;
  err += nc_def_var (dataFile, "g_psi_x", NC_DOUBLE, 1, &Q_d, &Gpx);
  
  // P_psi_x
  int Ppx;
  err += nc_def_var (dataFile, "P_psi_x", NC_DOUBLE, 1, &Q_d, &Ppx);

  // r/a
  double* rra = new double[NPSI];
  for (int j = 0; j < NPSI; j++)
    rra [j] = rP[j] /ra;
  int r_y;
  err += nc_def_var (dataFile, "r", NC_DOUBLE, 1, &Q_d, &r_y);
    
  // theta
  int T_d;
  err += nc_def_dim (dataFile, "k", NTHETA, &T_d);
  double* tth = new double[NTHETA];
  for (int j = 0; j < NTHETA; j++)
    tth [j] = th[j] /M_PI;
  int th_y;
  err += nc_def_var (dataFile, "theta", NC_DOUBLE, 1, &T_d, &th_y);
 
  // RRst
  double* DATA5 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA5[j + i*NTHETA] = gsl_matrix_get (RRst, i, j);

  int rrst_d[2], rrst;
  rrst_d[0] = Q_d;
  rrst_d[1] = T_d;
  err += nc_def_var (dataFile, "RR_st", NC_DOUBLE, 2, rrst_d, &rrst);

  // ZZst
  double* DATA6 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA6[j + i*NTHETA] = gsl_matrix_get (ZZst, i, j);

  int zzst;
  err += nc_def_var (dataFile, "ZZ_st", NC_DOUBLE, 2, rrst_d, &zzst);

  // RRr
  double* DATA7 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA7[j + i*NTHETA] = gsl_matrix_get (RRr, i, j);

  int dRdr;
  err += nc_def_var (dataFile, "dRdr", NC_DOUBLE, 2, rrst_d, &dRdr);

  // RRth
  double* DATA8 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA8[j + i*NTHETA] = gsl_matrix_get (RRth, i, j);

  int dRdth;
  err += nc_def_var (dataFile, "dRdth", NC_DOUBLE, 2, rrst_d, &dRdth);

  // ZZr
  double* DATA9 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA9[j + i*NTHETA] = gsl_matrix_get (ZZr, i, j);

  int dZdr;
  err += nc_def_var (dataFile, "dZdr", NC_DOUBLE, 2, rrst_d, &dZdr);

  // ZZth
  double* DATA10 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA10[j + i*NTHETA] = gsl_matrix_get (ZZth, i, j);

  int dZdth;
  err += nc_def_var (dataFile, "dZdth", NC_DOUBLE, 2, rrst_d, &dZdth);

  // Jac
  double* DATA11 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA11[j + i*NTHETA] = gsl_matrix_get (Jac, i, j);

  int jac;
  err += nc_def_var (dataFile, "Jac", NC_DOUBLE, 2, rrst_d, &jac);

  // Jax
  double* DATA12 = new double[NPSI*NTHETA];
  for (int i = 0; i < NPSI; i++)
    for (int j = 0; j < NTHETA; j++)
      DATA12[j + i*NTHETA] = gsl_matrix_get (Jax, i, j);

  int jax;
  err += nc_def_var (dataFile, "Jax", NC_DOUBLE, 2, rrst_d, &jax);

  // Write data
  err += nc_put_var_double (dataFile, PN,       PsiN);
  err += nc_put_var_double (dataFile, Q,        QP);
  err += nc_put_var_double (dataFile, Qx,       QX);
  err += nc_put_var_double (dataFile, G,        GP);
  err += nc_put_var_double (dataFile, F,        FP);
  err += nc_put_var_double (dataFile, Gp,       GPP);
  err += nc_put_var_double (dataFile, Px,       PP);
  err += nc_put_var_double (dataFile, Pp,       PPP);
  err += nc_put_var_double (dataFile, Gpx,      GPX);
  err += nc_put_var_double (dataFile, Ppx,      PPX);
  err += nc_put_var_double (dataFile, r_y,      rra);
  err += nc_put_var_double (dataFile, th_y,     tth);
  err += nc_put_var_double (dataFile, rrst,     DATA5);
  err += nc_put_var_double (dataFile, zzst,     DATA6);
  err += nc_put_var_double (dataFile, dRdr,     DATA7);
  err += nc_put_var_double (dataFile, dRdth,    DATA8);
  err += nc_put_var_double (dataFile, dZdr,     DATA9);
  err += nc_put_var_double (dataFile, dZdth,    DATA10);
  err += nc_put_var_double (dataFile, jac,      DATA11);
  err += nc_put_var_double (dataFile, jax,      DATA12);
  */

  /*
  delete[] rra;     delete[] tth;     
  delete[] DATA;    delete[] DATA1,   delete[] DATA2;   
  delete[] DATA5;   delete[] DATA6;   delete[] DATA7;
  delete[] DATA8;   delete[] DATA9;   delete[] DATA10;  delete[] DATA11;
  delete[] DATA12;
  */
 }

