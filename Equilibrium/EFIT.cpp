// EFIT.cpp

#include "Equilibrium.h"

// ###############################
// Function to calculate EFIT data
// ###############################
void Equilibrium::CalculateEFIT ()
{
  // ----------------------
  // Set physical constants
  // ----------------------
  double mu_0 = 4.*M_PI*1.e-7;

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  {
    string JSONFilename = "../Inputs/Equilibrium.json";
    json   JSONData     = ReadJSONFile (JSONFilename);

    EFIT  = JSONData["EFIT"] .get<int>    ();
    NRBOX = JSONData["NRBOX"].get<int>    ();
    NZBOX = JSONData["NZBOX"].get<int>    ();
    rc    = JSONData["rc"]   .get<double> ();
  }
  if (!EFIT)
    return;

  printf ("Calculating EFIT data:\n");
  printf ("NRBOX = %4d        NZBOX = %4d       rc = %10.3e\n", NRBOX, NZBOX, rc);

  {
    string JSONFilename = "../Inputs/TJ.json";
    json   JSONData     = ReadJSONFile (JSONFilename);
    
    B0EXP = JSONData["B0"].get<double> ();
    R0EXP = JSONData["R0"].get<double> ();
  }

  printf ("B0EXP = %10.3e  R0EXP = %10.3e\n", B0EXP, R0EXP);

  RAXIS   = R0EXP;
  ZAXIS   = 0.;
  CURRENT = epsa*epsa * (B0EXP * R0EXP /mu_0) * It[Nr];

  // ---------------
  // Allocate memory
  // ---------------
  PSI  = new double[NRBOX];
  PSIN = new double[NRBOX];
  rPSI = new double[NRBOX];
  T    = new double[NRBOX];
  TTp  = new double[NRBOX];
  P    = new double[NRBOX];
  Pp   = new double[NRBOX];
  Q    = new double[NRBOX];

  RBOUND   = new double[Nw+1];
  ZBOUND   = new double[Nw+1];
  RLIMITER = new double[5];
  ZLIMITER = new double[5];

  RGRID = new double[NRBOX];
  ZGRID = new double[NZBOX];
  PSIRZ = new double[NRBOX*NZBOX];
  rRZ   = new double[NRBOX*NZBOX];
  wRZ   = new double[NRBOX*NZBOX];
  cwRZ  = new double[NRBOX*NZBOX];
  swRZ  = new double[NRBOX*NZBOX];

  rPsispline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  PSIrspline = gsl_spline_alloc (gsl_interp_cspline, NRBOX);
  rPsiacc    = gsl_interp_accel_alloc ();
  PSIracc    = gsl_interp_accel_alloc ();

  // ------------------------
  // Set up Psi and PSI grids
  // ------------------------
  for (int i = 0; i <= Nr; i++)
    Psi[i] *= epsa*epsa * B0EXP * R0EXP*R0EXP;

  for (int i = 0; i <= Nr; i++)
    Psi[i] -= Psi[Nr];

  PSIAXIS  = Psi[0];
  PSIBOUND = Psi[Nr];

  for (int i = 0; i < NRBOX; i++)
    PSI[i] = PSIAXIS * (double (NRBOX - 1 - i) /double (NRBOX - 1));

  for (int i = 0; i < NRBOX; i++)
    PSIN[i] = (PSI[i] - PSIAXIS) /(PSIBOUND - PSIAXIS);

  // ---------------------------
  // Interpolate r onto Psi grid
  // ---------------------------
  gsl_spline_init (rPsispline, Psi, rr, Nr+1);

  // ----------------
  // Set up rPSI grid
  // ----------------
  rPSI[0]       = 0.;
  rPSI[NRBOX-1] = 1.;
  for (int i = 1; i < NRBOX - 1; i++)
    rPSI[i] = gsl_spline_eval (rPsispline, PSI[i], rPsiacc);
  
  // ---------------------------
  // Interpolate PSI onto r grid
  // ---------------------------
  gsl_spline_init (PSIrspline, rPSI, PSI, NRBOX);

  // ---------------------------
  // Calculate profile functions
  // ---------------------------
  double f1c   = 1./qc;
  double f3c   = - f1c * (HHfunc(2, 0) * HHfunc(2, 0) + VVfunc(2, 0) * VVfunc(2, 0));
  double p2ppc = - 2. * pc * mu;
  double g2pc  = - 2. * (f1c*f1c + p2ppc/2.);

  T  [0] = B0EXP * R0EXP;
  TTp[0] = B0EXP                      * g2pc  /(f1c + eps*eps * f3c);
  P  [0] = (B0EXP*B0EXP /mu_0)        * epsa*epsa * pc;
  Pp [0] = (B0EXP /mu_0 /R0EXP/R0EXP) * p2ppc /(f1c + eps*eps * f3c);
  Q  [0] = q2[0];
  
  for (int i = 1; i < NRBOX; i++)
    {
      double rp = rPSI[i];

      double f1  = Getf1 (rp);
      double f1p = Getf1p(rp);
      double p2  = Getp2 (rp);
      double p2p = Getp2p(rp);

      double g2  = gsl_spline_eval (g2spline, rp, g2acc);
      double f   = gsl_spline_eval (fspline,  rp, facc);
      double q   = gsl_spline_eval (q2spline, rp, q2acc);

      double g   = 1. + epsa*epsa * g2;
      double g2p = - p2p - f1*f1p/rp/rp;

      T  [i] = B0EXP * R0EXP              * g;
      TTp[i] = B0EXP                      * rp * g2p * g /f;
      P  [i] = (B0EXP*B0EXP /mu_0)        * epsa*epsa * p2;
      Pp [i] = (B0EXP /mu_0 /R0EXP/R0EXP) * rp * p2p     /f;
      Q  [i] = q;                             
    }

  // -------------------------
  // Calculate boundary points
  // -------------------------
  NPBOUND = Nw + 1;
  for (int i = 0; i < NPBOUND; i++)
    {
      RBOUND[i] = R0EXP * Rbound[i];
      ZBOUND[i] = R0EXP * Zbound[i];
    }

  // -----------------------------------------
  // Calculate bounding box and limiter points
  // -----------------------------------------
  double Rmin = 1.e6, Rmax = -1.e6, Zmin = 1.e6, Zmax = -1.e6;

  for (int i = 0; i < NPBOUND; i++)
    {
      if (RBOUND[i] < Rmin)
	Rmin = RBOUND[i];
      if (RBOUND[i] > Rmax)
	Rmax = RBOUND[i];
      if (ZBOUND[i] < Zmin)
	Zmin = ZBOUND[i];
      if (ZBOUND[i] > Zmax)
	Zmax = ZBOUND[i];
    }
  
  Rmin = Rmin - 0.05 * (Rmax - Rmin);
  Rmax = Rmax + 0.05 * (Rmax - Rmin);
  Zmin = Zmin - 0.05 * (Zmax - Zmin);
  Zmax = Zmax + 0.05 * (Zmax - Zmin);

  RBOXLFT = Rmin;
  RBOXLEN = Rmax - Rmin;
  ZOFF    = (Zmin + Zmax) /2.;
  ZBOXLEN = Zmax - Zmin;

  NLIMITER = 5;

  RLIMITER[0] = Rmin; ZLIMITER[0] = Zmin;
  RLIMITER[1] = Rmax; ZLIMITER[1] = Zmin;
  RLIMITER[2] = Rmax; ZLIMITER[2] = Zmax;
  RLIMITER[3] = Rmin; ZLIMITER[3] = Zmax;
  RLIMITER[4] = Rmin; ZLIMITER[4] = Zmin;

  // ------------------------------------
  // Set up normalized R and Z gridpoints
  // ------------------------------------
  for (int i = 0; i < NRBOX; i++)
    RGRID[i] = (RBOXLFT           + RBOXLEN * double (i) /double (NRBOX - 1)) /R0EXP;
  for (int j = 0; j < NZBOX; j++)
    ZGRID[j] = (ZOFF - ZBOXLEN/2. + ZBOXLEN * double (j) /double (NZBOX - 1)) /R0EXP;

  // -----------------------
  // Calculate PSIRZ on grid
  // -----------------------
  int cnt = 0;
  for (int i = 0; i < NRBOX; i++)
    {
      for (int j = 0; j < NZBOX; j++)
	{
	  double R = RGRID[i];
	  double Z = ZGRID[j];
	  
	  double r = sqrt ((R - 1.) * (R - 1.) + Z*Z) /epsa;
	  double w = atan2 (Z, 1. - R);
	  
	  for (int k = 0; k < 10; k++)
	    {
	      double RR = R - Getf_R (r, w);
	      double ZZ = Z - Getf_Z (r, w);
	      
	      r = sqrt ((RR - 1.) * (RR - 1.) + ZZ*ZZ) /epsa;
	      w = atan2 (ZZ, 1. - RR);
	    }
	  rRZ [cnt] = r;
	  wRZ [cnt] = w;
	  cwRZ[cnt] = cos (w);
	  swRZ[cnt] = sin (w);
	  
	  if (r > rc)
	    PSIRZ[cnt] = epsa*epsa * B0EXP*R0EXP*R0EXP * GetPSIvac (rc) * r*r/rc/rc;
	  else if (r >= 1.)
	    PSIRZ[cnt] = epsa*epsa * B0EXP*R0EXP*R0EXP * GetPSIvac (r);
	  else
	    PSIRZ[cnt] = gsl_spline_eval (PSIrspline, r, PSIracc);
	  
	  cnt++;
	}
      if (i%128 == 0)
	{
	  printf ("%04d ", i); fflush (stdout);
	}
    }
  printf ("\n");

  // --------------------------------------
  // Set up unnormalized R and Z gridpoints
  // --------------------------------------
  for (int i = 0; i < NRBOX; i++)
    RGRID[i] *= R0EXP;
  for (int j = 0; j < NZBOX; j++)
    ZGRID[j] *= R0EXP;

  // -------------------------------
  // Output EFIT data to netcdf file
  // -------------------------------
  printf ("Writing EFIT data to netcdf file Outputs/Equilibrium/EFIT.nc:\n");
  
  int pint[4];

  pint[0] = NRBOX;
  pint[1] = NZBOX;
  pint[2] = NPBOUND;
  pint[3] = NLIMITER;

  double preal[15];

  preal[0]  = RBOXLEN;
  preal[1]  = ZBOXLEN;
  preal[2]  = RBOXLFT;
  preal[3]  = ZOFF;
  preal[4]  = R0EXP;
  preal[5]  = B0EXP;
  preal[6]  = RAXIS;
  preal[7]  = ZAXIS;
  preal[8]  = PSIAXIS;
  preal[9]  = PSIBOUND;
  preal[10] = CURRENT;
  preal[11] = Rmin;
  preal[12] = Rmax;
  preal[13] = Zmin;
  preal[14] = Zmax;
  
  try
    {
      NcFile dataFile ("../Outputs/Equilibrium/EFIT.nc", NcFile::replace);

      NcDim i_d = dataFile.addDim ("Ni", 4);
      NcDim r_d = dataFile.addDim ("Nr", 15);
      NcDim p_d = dataFile.addDim ("Np", NRBOX);
      NcDim z_d = dataFile.addDim ("Nz", NZBOX);
      NcDim b_d = dataFile.addDim ("Nb", NPBOUND);
      NcDim l_d = dataFile.addDim ("Nl", NLIMITER);

      vector<NcDim> psi_d;
      psi_d.push_back (p_d);
      psi_d.push_back (z_d);
 
      NcVar i_x   = dataFile.addVar ("IntegerParameters", ncInt,    i_d);
      i_x.putVar (pint);
      NcVar r_x   = dataFile.addVar ("RealParameters",    ncDouble, r_d);
      r_x.putVar (preal); 
      NcVar pn_x  = dataFile.addVar ("PSI_N",             ncDouble, p_d);
      pn_x.putVar (PSIN);
      NcVar T_x   = dataFile.addVar ("T",                 ncDouble, p_d);
      T_x.putVar (T);
      NcVar P_x   = dataFile.addVar ("P",                 ncDouble, p_d);
      P_x.putVar (P);
      NcVar TTp_x = dataFile.addVar ("TTp",               ncDouble, p_d);
      TTp_x.putVar (TTp);
      NcVar Pp_x  = dataFile.addVar ("Pp",                ncDouble, p_d);
      Pp_x.putVar (Pp);
      NcVar PSI_x = dataFile.addVar ("PSI",               ncDouble, psi_d);
      PSI_x.putVar (PSIRZ);
      NcVar rx_x = dataFile.addVar  ("r",                 ncDouble, psi_d);
      rx_x.putVar (rRZ); 
      NcVar w_x = dataFile.addVar   ("w",                 ncDouble, psi_d);
      w_x.putVar (wRZ);
      NcVar cw_x = dataFile.addVar  ("cosw",              ncDouble, psi_d);
      cw_x.putVar (cwRZ);
      NcVar sw_x = dataFile.addVar  ("sinw",              ncDouble, psi_d);
      sw_x.putVar (swRZ);
      NcVar Q_x   = dataFile.addVar ("Q",                 ncDouble, p_d);
      Q_x.putVar (Q);
      NcVar R_x   = dataFile.addVar ("RBOUND",            ncDouble, b_d);
      R_x.putVar (RBOUND);
      NcVar Z_x   = dataFile.addVar ("ZBOUND",            ncDouble, b_d);
      Z_x.putVar (ZBOUND);
      NcVar rr_x  = dataFile.addVar ("RLIMITER",          ncDouble, l_d);
      rr_x.putVar (RLIMITER);
      NcVar zz_x  = dataFile.addVar ("ZLIMITER",          ncDouble, l_d);
      zz_x.putVar (ZLIMITER);
      NcVar rg_x  = dataFile.addVar ("RGRID",             ncDouble, p_d);
      rg_x.putVar (RGRID);
      NcVar zg_x  = dataFile.addVar ("ZGRID",             ncDouble, z_d);
      zg_x.putVar (ZGRID);
    }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file Outputs/Equilibrium/EFIT.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  // -------------------
  // Output PSI in ascii
  // -------------------
  FILE* file = OpenFilew ("../Outputs/Equilibrium/PsiSequential.txt");
  cnt = 0;
  for (int i = 0; i < NRBOX; i++)
    for (int j = 0; j < NRBOX; j++)
      {
	fprintf (file, "%4d %4d %17.9e\n", i, j, PSIRZ[cnt]);
	cnt++;
      }
  fclose (file);

  // --------
  // Clean up
  // --------
  delete[] PSI;    delete[] rPSI;   delete[] T;        delete[] TTp;      delete[] Pp;   
  delete[] RBOUND; delete[] ZBOUND; delete[] RLIMITER; delete[] ZLIMITER; delete[] PSIRZ;
  delete[] RGRID;  delete[] ZGRID;  delete[] P;        delete[] PSIN;     delete[] rRZ;
  delete[] wRZ;    delete[] Q;      delete[] cwRZ;     delete[] swRZ; 

  gsl_spline_free (rPsispline);    gsl_spline_free (PSIrspline);
  gsl_interp_accel_free (rPsiacc); gsl_interp_accel_free (PSIracc);
}
