// Rational.cpp

#include "TJ.h"

// #############################################
// Target function for finding rational surfaces
// #############################################
double TJ::RootFindF (double r)
{
  double q = Getq (r);
  
  return qval - q;
}

// ###################################
// Function to find rational surfaces.
// Assumes monotonic q-profile.
// ###################################
void TJ::FindRational ()
{
  // .....................................
  // Determine number of rational surfaces
  // .....................................
  double q0   = Getq (EPS);
  double qa   = Getq (1.);
  int    mmin = int (ntor * q0) + 1;
  int    mmax = int (ntor * qa);
  nres        = mmax - mmin + 1;

  // Abort calculation if no rational surfaces in plasma
  if (nres <= 0)
    {
      printf ("TJ::FindRational: Error - no rational surfaces in plasma\n");
      exit (1);
    }

  // .....................................
  // Calculate rational surface quantities
  // .....................................
  mres   = new int   [nres];
  qres   = new double[nres];
  rres   = new double[nres];
  Pres   = new double[nres];
  qerr   = new double[nres];
  sres   = new double[nres];
  gres   = new double[nres];
  DIres  = new double[nres];
  DRres  = new double[nres];
  nuLres = new double[nres];
  nuSres = new double[nres];
  Jres   = new int   [nres];
  Flarge = new double[nres];
  Fsmall = new double[nres];
  neres  = new double[nres];
  nepres = new double[nres];
  Teres  = new double[nres];
  Tepres = new double[nres];
 
  for (int i = 0; i < nres; i++)
    {
      mres[i]  = mmin + i;
      qres[i]  = double (mres[i]) /ntor;

      qval     = qres[i];
      rres[i]  = RootFind (EPS, 1.);
      
      qerr  [i] = fabs (Getq (rres[i]) - qres[i]);
      sres  [i] = Gets (rres[i]);
      gres  [i] = 1. + epsa*epsa * Getg2 (rres[i]);
      DIres [i] = GetDI (rres[i]);
      DRres [i] = GetDR (rres[i]);
      Pres  [i] = GetPsiN (rres[i]);
      neres [i] = Getne (rres[i]);
      nepres[i] = Getnep (rres[i]);
      Teres [i] = GetTe (rres[i]);
      Tepres[i] = GetTep (rres[i]);

      if (DIres[i] > 0.)
	{
	  printf ("TJ::FindRational - rational surface %3d unstable to ideal interchange modes\n", i);
	  exit (1);
	}
      
      nuLres[i] = 0.5 - sqrt (- DIres[i]);
      nuSres[i] = 0.5 + sqrt (- DIres[i]);

      Flarge[i] = GetFlarge (rres[i], mres[i]);
      Fsmall[i] = GetFsmall (rres[i], mres[i]);
    }

  for (int k = 0; k < nres; k++)
    for (int j = 0; j < J; j++)
      if (MPOL[j] == mres[k]) 
	Jres[k] = j;

  printf ("Rational surface data:\n");
  for (int i = 0; i < nres; i++)
    {
    printf ("m = %3d PsiN = %10.3e, r   = %10.3e s = %10.3e DI = %10.3e DR = %10.3e nuL = %10.3e nuS = %10.3e |q-qs| = %10.3e\n",
	    mres[i], Pres[i], rres[i], sres[i], DIres[i], DRres[i], nuLres[i], nuSres[i], qerr[i]);
    if (TEMP)
      printf ("        Tep  = %10.3e  nep = %10.3e\n",
	      Tepres[i], nepres[i]);
    }
  
  // Abort calculation if resonant mode numbers outside range of included poloidal harmonics
  if (mres[0] < MMIN || mres[nres-1] > MMAX)
    {
      printf ("TJ::FindRational: Error - resonant mode numbers do not lie in range MMIN to MMAX\n");
      exit (1);
    }
}

// ###################################
// Function to get resonant layer data
// ###################################
void TJ::GetLayerData ()
{
  // ...............
  // Allocate memory
  // ...............
  S13res = new double[nres];
  taures = new double[nres];
  ieres  = new double[nres];
  QEres  = new double[nres];
  Qeres  = new double[nres];
  Qires  = new double[nres];
  Dres   = new double[nres];
  Pmres  = new double[nres];
  Peres  = new double[nres];
  Dcres  = new double[nres];

  // .........................
  // Define physical constants
  // .........................
  double e    = 1.602176634e-19;
  double mp   = 1.672621925956e-27;
  double me   = 9.1093837139e-31;
  double mu0  = 4.*M_PI*1.e-7;
  double eps0 = 8.8541878188e-12;

  // ..........................
  // Calculate layer quantities
  // ..........................
  for (int k = 0; k < nres; k++)
    {
      double g2k = gsl_spline_eval (g2spline, rres[k], g2acc);
      double p2k = gsl_spline_eval (p2spline, rres[k], p2acc);
      double ppk = gsl_spline_eval (ppspline, rres[k], ppacc);

      double nek  = neres[k];
      double Tek  = Teres[k];
      double Lnk  = 24. + 3.*log(10.) - 0.5 * log(nek) + log(Tek);
      double teek = 6.*sqrt(2.)*pow(M_PI, 1.5) * eps0*eps0 * sqrt(me) * pow(Tek, 1.5)
	/Lnk /pow(e, 2.5) /nek;
      double sigk = ((sqrt(2.) + 13.*Zeff/4.)/Zeff /(sqrt(2.) + Zeff)) * nek *e*e * teek /me;
      double gk   = 1. + epsa*epsa * g2k;
      double Lsk  = R0 * qres[k] /sres[k];
      double VAk  = B0 * gk /sqrt (mu0 * nek * Mion * mp);
      double dik  = sqrt (Mion * mp /nek /e/e /mu0);
      double bek  = (5./3.) * epsa*epsa * p2k /gk/gk;
      double dbk  = sqrt (bek /(1. + bek)) * dik /apol /rres[k];
      double wak  = double(mres[k]) * B0 * ppk /mu0 /R0/R0 /e /nek /gk /rres[k];
      double tHk  = Lsk /double(mres[k]) /VAk;
      double tRk  = mu0 * apol*apol * rres[k]*rres[k] * sigk;
      double tPk  = apol*apol * rres[k]*rres[k] /Chip;

      double vTe  = sqrt (2. * Tek*e /me);
      double Chps = 1.581 * teek * vTe*vTe /(1. + 0.2535 * Zeff);
      double Wd   = 0.1;

      for (int i = 0; i < ITERMAX; i++)
	{
	  double Chpl = 2.*R0 * vTe /sqrt(M_PI) /double(NTOR) /sres[k] /Wd;
	  double Chpa = Chps * Chpl /(Chps + Chpl);
	  Wd          = sqrt(8.) * pow (Chip /Chpa, 0.25) /sqrt (epsa * rres[k] * sres[k] * double(NTOR));
	}

      S13res[k] = pow (tRk/tHk, 1./3.);
      taures[k] = S13res[k] * tHk;
      ieres [k] = 0.5;
      Qeres [k] = - ieres[k]      * taures[k] * wak;
      Qires [k] = (1. - ieres[k]) * taures[k] * wak;
      QEres [k] = 0.;
      Dres  [k] = S13res[k] * sqrt (ieres[k]) * dbk;
      Pmres [k] = tRk /tPk;
      Peres [k] = tRk /tPk;
      Dcres [k] = - sqrt(2.) * pow(M_PI, 1.5) * DRres[k] /Wd;
    }

  // .......................
  // Output layer quantities
  // .......................
  printf ("Resonant layer data:\n");
  for (int k = 0; k < nres; k++)
    printf ("m = %3d r = %9.2e Te = %9.2e S13 = %9.2e tau = %9.2e Qe = %9.2e D = %9.2e P = %9.2e Delta_c = %9.2e\n",
	    mres[k], rres[k], Teres[k], S13res[k], taures[k], Qeres[k], Dres[k], Pmres[k], Dcres[k]);
}
