// ReadEquilibrium.cpp

#include "TJ.h"

// #############################################################
// Function to read equilibrium data from Outputs/Equilibrium.nc
// #############################################################
void TJ::ReadEquilibrium ()
{
  try
    {
      // .........................................
      // Read in equilibrium data from NETCDF file
      // .........................................
      NcFile dataFile ("Plots/Equilibrium.nc", NcFile::read);

      NcVar p_x = dataFile.getVar ("para");
      NcDim p_d = p_x.getDim (0);

      int     Np   = p_d.getSize ();
      double* para = new double[Np];

      p_x.getVar(para);
      epsa = para[0];
      sa   = para[1];
      G1   = para[2];
      G2   = para[3];
      
      NcVar r_x   = dataFile.getVar ("r");
      NcVar pp_x  = dataFile.getVar ("pp");
      NcVar ppp_x = dataFile.getVar ("ppp");
      NcVar q_x   = dataFile.getVar ("q");
      NcVar s_x   = dataFile.getVar ("s");
      NcVar s2_x  = dataFile.getVar ("s2");
      NcVar S1_x  = dataFile.getVar ("S1");
      NcVar P1_x  = dataFile.getVar ("P1");
      NcVar P2_x  = dataFile.getVar ("P2");
      NcVar P3_x  = dataFile.getVar ("P3");
      NcDim r_d   = r_x.getDim (0);

      Nr  = r_d.getSize () - 1;
      rr  = new double[Nr+1];
      pp  = new double[Nr+1];
      ppp = new double[Nr+1];
      q   = new double[Nr+1];
      s   = new double[Nr+1];
      s2  = new double[Nr+1];
      S1  = new double[Nr+1];
      P1  = new double[Nr+1];
      P2  = new double[Nr+1];
      P3  = new double[Nr+1];

      r_x.  getVar (rr);
      pp_x. getVar (pp);
      ppp_x.getVar (ppp);
      q_x  .getVar (q);
      s_x  .getVar (s);
      s2_x .getVar (s2);
      S1_x .getVar (S1);
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
      
      delete[] para;   delete[] Hndata; delete[] Hnpdata; delete[] Vndata; delete[] Vnpdata;
      delete[] RRdata; delete[] ZZdata; delete[] rrdata;  delete[] ttdata;
    }
  catch (NcException& e)
     {
       e.what ();
       printf ("Error reading Inputs/Equilbrium.nc\n");
       exit (1);
     }

  // .....................................................
  // Allocate memory for interpolation of equilibrium data
  // .....................................................
  ppspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  pppspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  qspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  sspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  s2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  S1spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P1spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P3spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

  ppacc     = gsl_interp_accel_alloc ();
  pppacc    = gsl_interp_accel_alloc ();
  qacc      = gsl_interp_accel_alloc ();
  sacc      = gsl_interp_accel_alloc ();
  s2acc     = gsl_interp_accel_alloc ();
  S1acc     = gsl_interp_accel_alloc ();
  P1acc     = gsl_interp_accel_alloc ();
  P2acc     = gsl_interp_accel_alloc ();
  P3acc     = gsl_interp_accel_alloc ();

  HHspline  = new gsl_spline* [Ns+1];
  VVspline  = new gsl_spline* [Ns+1];
  HPspline  = new gsl_spline* [Ns+1];
  VPspline  = new gsl_spline* [Ns+1];

  HHacc     = new gsl_interp_accel* [Ns+1];
  VVacc     = new gsl_interp_accel* [Ns+1];
  HPacc     = new gsl_interp_accel* [Ns+1];
  VPacc     = new gsl_interp_accel* [Ns+1];

  for (int n = 0; n <= Ns; n++)
    {
      HHspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VVspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      HPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
      VPspline[n] = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

      HHacc[n]    = gsl_interp_accel_alloc ();
      VVacc[n]    = gsl_interp_accel_alloc ();
      HPacc[n]    = gsl_interp_accel_alloc ();
      VPacc[n]    = gsl_interp_accel_alloc ();
    }

  // ............................
  // Interpolate equilibrium data
  // ............................
  gsl_spline_init (ppspline,  rr, pp,  Nr+1);
  gsl_spline_init (pppspline, rr, ppp, Nr+1);
  gsl_spline_init (qspline,   rr, q,   Nr+1);
  gsl_spline_init (sspline,   rr, s,   Nr+1);
  gsl_spline_init (s2spline,  rr, s2,  Nr+1);
  gsl_spline_init (S1spline,  rr, S1,  Nr+1);
  gsl_spline_init (P1spline,  rr, P1,  Nr+1);
  gsl_spline_init (P2spline,  rr, P2,  Nr+1);
  gsl_spline_init (P3spline,  rr, P3,  Nr+1);

  double* data = new double[Nr+1];

  for (int i = 0; i <= Nr; i++)
    data[i] = HHfunc(1, i);
  gsl_spline_init (HHspline[1], rr, data, Nr+1);

  for (int i = 0; i <= Nr; i++)
    data[i] = HPfunc(1, i);
  gsl_spline_init (HPspline[1], rr, data, Nr+1);

  for (int n = 2; n <= Ns; n++)
    {
      for (int i = 0; i <= Nr; i++)
	data[i] = HHfunc(n, i);
      gsl_spline_init (HHspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = HPfunc(n, i);
      gsl_spline_init (HPspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VVfunc(n, i);
      gsl_spline_init (VVspline[n], rr, data, Nr+1);

      for (int i = 0; i <= Nr; i++)
	data[i] = VPfunc(n, i);
      gsl_spline_init (VPspline[n], rr, data, Nr+1);
    }
  delete[] data;

  // Output equilibrium data
  printf ("Equilibrium data:\n");
  printf ("epsa = %10.3e q0 = %10.3e qa = %10.3e sa = %10.3e G1 = %10.3e G2 = %10.3e\n",
	  epsa, Getq (0.), Getq (1.), sa, G1, G2);
  printf ("n = %3d Hna = %10.3e Vna = %10.3e\n", 1, GetHn (1, 1.), 0.);
  for (int n = 2; n <= Ns; n++)
    if (GetHn (n, 1.) > 1.e-15 || GetVn (n, 1.) > 1.e-15)
      printf ("n = %3d Hna = %10.3e Vna = %10.3e\n", n, GetHn (n, 1.), GetVn (n, 1.));
}

// ####################################################
// Function to calculate metric data on plasma boundary
// ####################################################
void TJ::CalculateMetric ()
{
  // ...............
  // Allocate memory
  // ...............
  cmu    = new double[Nw+1];
  eeta   = new double[Nw+1];
  ceta   = new double[Nw+1];
  seta   = new double[Nw+1];
  R2grgz = new double[Nw+1];
  R2grge = new double[Nw+1];

  Rbspline = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Zbspline = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Rbacc    = gsl_interp_accel_alloc ();
  Zbacc    = gsl_interp_accel_alloc ();
 
  Rrzspline = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Rrespline = gsl_spline_alloc (gsl_interp_cspline_periodic, Nw+1);
  Rrzacc    = gsl_interp_accel_alloc ();
  Rreacc    = gsl_interp_accel_alloc ();
 
  // .....................
  // Calculate metric data
  // .....................
  for (int i = 0; i <= Nw; i++)
    {
      double R  = Rbound[i];
      double Z  = Zbound[i];
      double Rt = dRdthe[i];
      double Zt = dZdthe[i];

      double z   = GetCoshMu (R, Z);
      double et  = GetEta    (R, Z);
      double cet = cos (et);
      double set = sin (et);

      double muR = 1. - z * cet;
      double muZ = - sqrt (z*z - 1.) * set;
      double etR = - sqrt (z*z - 1.) * set;
      double etZ = z * cet - 1.;

      cmu [i] = z;
      eeta[i] = et /M_PI;
      ceta[i] = cet;
      seta[i] = set;
    
      R2grgz[i] = R * sqrt (z*z - 1.) * (Rt * muZ - Zt * muR);
      R2grge[i] = R                   * (Rt * etZ - Zt * etR);
    }

  // .......................
  // Interpolate metric data
  // .......................
  gsl_spline_init (Rbspline,  tbound, Rbound, Nw+1);
  gsl_spline_init (Zbspline,  tbound, Zbound, Nw+1);
  gsl_spline_init (Rrzspline, tbound, R2grgz, Nw+1);
  gsl_spline_init (Rrespline, tbound, R2grge, Nw+1);
}
