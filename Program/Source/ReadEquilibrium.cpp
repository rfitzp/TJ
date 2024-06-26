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
      Nw = w_d.getSize ();

      RR.resize     (Nf, Nw);
      ZZ.resize     (Nf, Nw);
      rvals.resize  (Nf, Nw);
      thvals.resize (Nf, Nw);

      double* RRdata = new double[Nf*Nw];
      double* ZZdata = new double[Nf*Nw];
      double* rrdata = new double[Nf*Nw];
      double* ttdata = new double[Nf*Nw];

      RR_x.getVar (RRdata);
      ZZ_x.getVar (ZZdata);
      rr_x.getVar (rrdata);
      tt_x.getVar (ttdata);

      rf = new double[Nf];
      for (int n = 0; n < Nf; n++)
	rf[n] = rrdata[n*Nw];
      
      for (int n = 0; n < Nf; n++)
	for (int i = 0; i < Nw; i++)
	  {
	    RR    (n, i) = RRdata[i + n*Nw];
	    ZZ    (n, i) = ZZdata[i + n*Nw];
	    rvals (n, i) = rrdata[i + n*Nw];
	    thvals(n, i) = ttdata[i + n*Nw];
	  }
      
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
  printf ("epsa = %10.3e q0 = %10.3e qa = %10.3e H1a = %10.3e sa = %10.3e G1 = %10.3e G2 = %10.3e\n",
	  epsa, Getq (0.), Getq (1.), GetHn (1, 1.), sa, G1, G2);
  for (int n = 2; n <= Ns; n++)
    printf ("n = %3d Hna = %10.3e Vna = %10.3e\n", n, GetHn (n, 1.), GetVn (n, 1.));
}
