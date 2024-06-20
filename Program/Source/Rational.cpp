// Rational.cpp

#include "TJ.h"

// #############################################
// Target function for finding rational surfaces
// #############################################
double TJ::Feval (double r)
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
  qerr   = new double[nres];
  sres   = new double[nres];
  DIres  = new double[nres];
  nuLres = new double[nres];
  nuSres = new double[nres];
  Jres   = new int   [nres];
  
  for (int i = 0; i < nres; i++)
    {
      mres[i]  = mmin + i;
      qres[i]  = double (mres[i]) /ntor;

      qval     = qres[i];
      rres[i]  = RootFind ();
      
      qerr [i]  = fabs (Getq (rres[i]) - qres[i]);
      sres [i]  = Gets (rres[i]);
      DIres[i]  = GetDI (rres[i]);

      if (DIres[i] > 0.)
	{
	  printf ("TJ::FindRational - rational surface %3d unstable to ideal interchange modes\n", i);
	  exit (1);
	}
      
      nuLres[i] = 0.5 - sqrt(-DIres[i]);
      nuSres[i] = 0.5 + sqrt(-DIres[i]);
    }

  for (int k = 0; k < nres; k++)
    for (int j = 0; j < J; j++)
      if (MPOL[j] == mres[k]) 
	Jres[k] = j;

  printf ("Rational surface data:\n");
  for (int i = 0; i < nres; i++)
    printf ("mpol = %3d rs = %10.3e s = %10.3e DI = %10.3e nuL = %10.3e nuS = %10.3e res(q-qres) = %10.3e\n",
	    mres[i], rres[i], sres[i], DIres[i], nuLres[i], nuSres[i], qerr[i]);

  // Abort calculation if resonant mode numbers outside range of included poloidal harmonics
  if (mres[0] < MMIN || mres[nres-1] > MMAX)
    {
      printf ("TJ::FindRational: Error - resonant mode numbers do not lie in range MMIN to MMAX\n");
      exit (1);
    }
}
