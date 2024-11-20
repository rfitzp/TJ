// ReadEquilibrium.cpp

#include "TJ.h"

// ###########################################################
// Function to read equilibrium data from Plots/Equilibrium.nc
// ###########################################################
void TJ::ReadEquilibrium ()
{
  // ......................................
  // Read equilibrium data from netcdf file
  // ......................................
  ReadNetcdf ();
  
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

// ###############################
// Function to calculate wall data
// ###############################
void TJ::CalcWall ()
{
  // ...............
  // Allocate memory
  // ...............
  Rwm = new double[Nw+1];
  Zwm = new double[Nw+1];
  Rwp = new double[Nw+1];
  Zwp = new double[Nw+1];

  Ipw.resize   (J, J);
  Iqw.resize   (J, J);
  Jpw.resize   (J, J);
  Jqw.resize   (J, J);
  kw.resize    (J, J);
  Ihpw.resize  (J, J);
  Ihqw.resize  (J, J);
  Iw.resize    (J, J);
  Iwher.resize (J, J);
  Iwant.resize (J, J);
  Jqpw.resize  (J, J);
  Jw.resize    (J, J);
  Kw.resize    (J, J);
  Kwher.resize (J, J);
  Kwant.resize (J, J);

  // ................................
  // Calculate zwm, zwp, swm, and swp
  // ................................
  zwm = 1. /bw /epsa;
  swm = sqrt (zwm*zwm - 1.);

  Dw  = epsa * dw * bw;
  zwp = cosh (acosh(zwm) - swm * Dw);
  swp = sqrt (zwp*zwp - 1.);

  // ........................................................
  // Calculate coordinates of inner and outer wall boundaries
  // ........................................................
  for (int i = 0; i <= Nw; i++)
    {
      double eta = double (i) * 2.*M_PI /double (Nw);
      
      Rwm[i] = swm      /(zwm - cos(eta));
      Zwm[i] = sin(eta) /(zwm - cos(eta));
      Rwp[i] = swp      /(zwp - cos(eta));
      Zwp[i] = sin(eta) /(zwp - cos(eta));
    }

  // .......................
  // Calculate wall matrices
  // .......................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	int    M   = MPOL[j];
	int    MP  = MPOL[jp];
	int    MMP = abs (MP);
	    
	double m   = mpol[j];
	double mp  = mpol[jp];

	double Ptor  = NormToroidalP    (NTOR, MMP, zwm);
	double Ptorz = NormToroidaldPdz (NTOR, MMP, zwm);
	double Qtor  = NormToroidalQ    (NTOR, MMP, zwm);
	double Qtorz = NormToroidaldQdz (NTOR, MMP, zwm);

	if (M == MP)
	  {
	    Ipw(j, jp) = (0.5 * Ptor + zwm * Ptorz);
	    Iqw(j, jp) = (0.5 * Qtor + zwm * Qtorz);
	    Jpw(j, jp) = (0.625 + (0.5 + zwm*zwm) * (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Ptor;
	    Jqw(j, jp) = (0.625 + (0.5 + zwm*zwm) * (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Qtor;
	    kw (j, jp) = zwm;
	  }
	else if (M == MP + 1 || M == MP - 1)
	  {
	    Ipw(j, jp) = - 0.5 * Ptorz;
	    Iqw(j, jp) = - 0.5 * Qtorz;
	    Jpw(j, jp) = - zwm * (0.25 + (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Ptor;
	    Jqw(j, jp) = - zwm * (0.25 + (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Qtor;
	    kw (j, jp) = - 0.5;
	  }
	else if (M == MP + 2 || M == MP - 2)
	  {
	    Ipw(j, jp) = 0.;
	    Iqw(j, jp) = 0.;
	    Jpw(j, jp) = (- 0.0625 + 0.25 * (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Ptor;
	    Jqw(j, jp) = (- 0.0625 + 0.25 * (mp*mp + ntor*ntor /(zwm*zwm - 1.))) * Qtor;
	    kw (j, jp) = 0.;
	  }
	else
	  {
	    Ipw(j, jp) = 0.;
	    Iqw(j, jp) = 0.;
	    Jpw(j, jp) = 0.;
	    Jqw(j, jp) = 0.;
	    kw (j, jp) = 0.;
	  }
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Ihpw(j, jp) = 0.;
	Ihqw(j, jp) = 0.;

	for (int jpp = 0; jpp < J; jpp++)
	  {
	    Ihpw(j, jp) += kw(j, jpp) * Ipw(jpp, jp);
	    Ihqw(j, jp) += kw(j, jpp) * Iqw(jpp, jp);
	  }
      }

  // ...................
  // Calculate Iw-matrix
  // ...................
  SolveLinearSystem (Ihqw, Iw, Ihpw);
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Iw(j, jp) *= -1.;

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Iwher(j, jp) = 0.5 * (Iw(j, jp) + Iw(jp, j));
	Iwant(j, jp) = 0.5 * (Iw(j, jp) - Iw(jp, j));
      }

  // .............................
  // Calculate Iw-matrix residuals
  // .............................
  double Ihmax = 0., Iamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ihval = fabs (Iwher(j, jp));
	double iaval = fabs (Iwant(j, jp));

	if (ihval > Ihmax)
	  Ihmax = ihval;
	if (iaval > Iamax)
	  Iamax = iaval;	
      }

  // ...................
  // Calculate Jw-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Jqpw(j, jp) = Jpw(j, jp);

	for (int jpp = 0; jpp < J; jpp++)
	  Jqpw(j, jp) += Jqw(j, jpp) * Iw(jpp, jp);
      }

  SolveLinearSystem (Ihpw, Jw, Jqpw);

  // ...................
  // Calculate Kw-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Kw(j, jp) = 0.;

	for (int jpp = 0; jpp < J; jpp++)
	  Kw(j, jp) += Iw(j, jpp) * Jw(jpp, jp);
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Kwher(j, jp) = 0.5 * (Kw(j, jp) + Kw(jp, j));
	Kwant(j, jp) = 0.5 * (Kw(j, jp) - Kw(jp, j));
      }

  // .............................
  // Calculate Kw-matrix residuals
  // .............................
  double Khmax = 0., Kamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double khval = fabs (Kwher(j, jp));
	double kaval = fabs (Kwant(j, jp));

	if (khval > Khmax)
	  Khmax = khval;
	if (kaval > Kamax)
	  Kamax = kaval;	
      }

  printf ("Wall matrix Hermitian test residual: %10.4e %10.4e\n", Iamax/Ihmax, Kamax/Khmax);
}

// ##############################
// Function to read RMP coil data
// ##############################
void TJ::ReadCoils ()
{
  // ..........................................
  // Read shaping data from file Inputs/TJ.json
  // ..........................................
  string         JSONFilename = "../Inputs/TJ.json";
  json           JSONData     = ReadJSONFile (JSONFilename);
  vector<double> etacoil, icoil;

  for (const auto& number : JSONData["etacoil"])
    {
      etacoil.push_back (number.get<double> ());
    }
  for (const auto& number : JSONData["Icoil"])
    {
      icoil.push_back (number.get<double> ());
    }

  if (etacoil.size() != icoil.size())
    {
      printf ("TJ:: Error reading etacoil and Icoil arrays must be same size\n");
      exit (1);
    }
  ncoil = etacoil.size();

  Rcoil   = new double[ncoil];
  Zcoil   = new double[ncoil];
  Icoil   = new double[ncoil];
  
  for (int i = 0; i < ncoil; i++)
    {
      Rcoil[i] = swp                  /(zwp - cos(etacoil[i]*M_PI));
      Zcoil[i] = sin(etacoil[i]*M_PI) /(zwp - cos(etacoil[i]*M_PI));
      Icoil[i] = icoil[i];
     }

  printf ("RMP coil data:\n");
  for (int i = 0; i < ncoil; i++)
    printf ("eta_coil/M_PI = %10.3e Rcoil = %10.3e Zcoil = %10.3e Icoil = %10.3e\n", etacoil[i], Rcoil[i], Zcoil[i], Icoil[i]);
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
