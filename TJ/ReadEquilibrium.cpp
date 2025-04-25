// ReadEquilibrium.cpp

#include "TJ.h"

// #########################################################################
// Function to read equilibrium data from Outputs/Equilibrium/Equilibrium.nc
// #########################################################################
void TJ::ReadEquilibrium ()
{
  // ......................................
  // Read equilibrium data from netcdf file
  // ......................................
  ReadEquilibriumNetcdf ();
  
  // .....................................................
  // Allocate memory for interpolation of equilibrium data
  // .....................................................
  Pspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  fspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  g2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  p2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  ppspline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  pppspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  qspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  sspline   = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  s2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  s0spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  S1spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  S2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  S3spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P1spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P2spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  P3spline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  nespline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Tespline  = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  nepspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);
  Tepspline = gsl_spline_alloc (gsl_interp_cspline, Nr+1);

  Pacc      = gsl_interp_accel_alloc ();
  facc      = gsl_interp_accel_alloc ();
  g2acc     = gsl_interp_accel_alloc ();
  p2acc     = gsl_interp_accel_alloc ();
  ppacc     = gsl_interp_accel_alloc ();
  pppacc    = gsl_interp_accel_alloc ();
  qacc      = gsl_interp_accel_alloc ();
  sacc      = gsl_interp_accel_alloc ();
  s2acc     = gsl_interp_accel_alloc ();
  s0acc     = gsl_interp_accel_alloc ();
  S1acc     = gsl_interp_accel_alloc ();
  S2acc     = gsl_interp_accel_alloc ();
  S3acc     = gsl_interp_accel_alloc ();
  P1acc     = gsl_interp_accel_alloc ();
  P2acc     = gsl_interp_accel_alloc ();
  P3acc     = gsl_interp_accel_alloc ();
  neacc     = gsl_interp_accel_alloc ();
  Teacc     = gsl_interp_accel_alloc ();
  nepacc    = gsl_interp_accel_alloc ();
  Tepacc    = gsl_interp_accel_alloc ();

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
  gsl_spline_init (Pspline,   rr, PsiN, Nr+1);
  gsl_spline_init (fspline,   rr, f,    Nr+1);
  gsl_spline_init (g2spline,  rr, g2,   Nr+1);
  gsl_spline_init (p2spline,  rr, p2,   Nr+1);
  gsl_spline_init (ppspline,  rr, pp,   Nr+1);
  gsl_spline_init (pppspline, rr, ppp,  Nr+1);
  gsl_spline_init (qspline,   rr, q,    Nr+1);
  gsl_spline_init (sspline,   rr, s,    Nr+1);
  gsl_spline_init (s2spline,  rr, s2,   Nr+1);
  gsl_spline_init (s0spline,  rr, s0,   Nr+1);
  gsl_spline_init (S1spline,  rr, S1,   Nr+1);
  gsl_spline_init (S2spline,  rr, S2,   Nr+1);
  gsl_spline_init (S3spline,  rr, S3,   Nr+1);
  gsl_spline_init (P1spline,  rr, P1,   Nr+1);
  gsl_spline_init (P2spline,  rr, P2,   Nr+1);
  gsl_spline_init (P3spline,  rr, P3,   Nr+1);
  gsl_spline_init (nespline,  rr, ne,   Nr+1);
  gsl_spline_init (Tespline,  rr, Te,   Nr+1);
  gsl_spline_init (nepspline, rr, nep,  Nr+1);
  gsl_spline_init (Tepspline, rr, Tep,  Nr+1);

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

  // .......................
  // Output equilibrium data
  // .......................
  apol = epsa * R0;
  printf ("Plasma equilibrium data:\n");
  printf ("epsa = %10.3e q0 = %10.3e qa = %10.3e sa = %10.3e apol = %10.3e\n",
	  epsa, Getq (0.), Getq (1.), sa, apol);
  printf ("n = %3d Hna = %10.3e Vna = %10.3e\n", 1, GetHn (1, 1.), 0.);
  for (int n = 2; n <= Ns; n++)
    if (fabs(GetHn (n, 1.)) > 1.e-15 || fabs(GetVn (n, 1.)) > 1.e-15)
      printf ("n = %3d Hna = %10.3e Vna = %10.3e\n", n, GetHn (n, 1.), GetVn (n, 1.));
}

// ####################################################
// Function to calculate metric data on plasma boundary
// ####################################################
void TJ::CalculateMetricBoundary ()
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

// ##############################################################
// Function to read island data from Outputs/Island/Islandxxxx.nc
// ##############################################################
void TJ::ReadIsland ()
{
  // ...............
  // Allocate memory
  // ...............
  deltaTh.resize (nres, Nh, NX);
  T0inf     = new double           [nres];
  XX        = new double           [NX];
  dThspline = new gsl_spline*      [nres*Nh];
  dThacc    = new gsl_interp_accel*[nres*Nh];

  for (int k = 0; k < nres; k++)
    for (int n = 0; n < Nh; n++)
      {
	dThspline[k*Nh + n] = gsl_spline_alloc (gsl_interp_cspline, NX);
	dThacc   [k*Nh + n] = gsl_interp_accel_alloc ();
      }

  // ....................................
  // Run Island for each rational surface
  // ....................................
  for (int k = 0; k < nres; k++)
    {
	Island island (k, delta[k]);

	island.Solve (0);
    }

  // ...........................................................................
  // Read island harmonic data for each rational surface from Island netcdf file
  // ...........................................................................
  for (int k = 0; k < nres; k++)
    ReadIslandNetcdf (k);
  
  // ................................
  // Interpolate island harmonic data
  // ................................
  double* data = new double[NX];
  
  for (int k = 0; k < nres; k++)
    {
      for (int n = 0; n < Nh; n++)
	{
 	  for (int i = 0; i < NX; i++)
	    data[i] = deltaTh (k, n, i);

	  gsl_spline_init (dThspline[k*Nh + n], XX, data, NX);
	}
    }

  delete[] data;
}
  
