// RMP.cpp

#include "TJ.h"

// #########################################################
// Function to calculate resonant magnetic perturbation data
// #########################################################
void TJ::CalculateResonantMagneticPerturbation ()
{
  // ...............
  // Allocate memory
  // .............
  Upsilon = new complex<double>[J];
  Lambda  = new complex<double>[J];
  Chi     = new complex<double>[J];

  Psirmp.resize(J, NDIAG);
  Zrmp  .resize(J, NDIAG);

  Psirmps = new complex<double>[Nw+1];
  Psixs   = new complex<double>[Nw+1];

  // ........................
  // Calculate Upsilon-vector
  // ........................
  SolveLinearSystem (Xmat, Upsilon, Xi);

  // .......................
  // Calculate Lambda-vector
  // .......................
  for (int k = 0; k < nres; k++)
    {
      complex<double> sum = complex<double> (0., 0.);

      for (int j = 0; j < J; j++)
	sum += Pia(k, j) * Upsilon[j];

      Lambda[k] = sum;
    }

  // ....................
  // Calculate Chi-vector
  // ....................
  for (int k = 0; k < nres; k++)
    {
      complex<double> sum = complex<double> (0., 0.);

      for (int kp = 0; kp < nres; kp++)
	sum += Emat(k, kp) * Lambda[kp];
      
      Chi[k] = sum;
    }

  // .................
  // Output Chi-vector
  // .................
  printf ("Chi vector:\n");
  for (int k = 0; k < nres; k++)
    printf ("Rational surface %2d: Chi = (%10.3e, %10.3e) |Chi| = %10.3e\n",
	    k+1, real(Chi[k]), imag(Chi[k]), abs(Chi[k]));

  // .........................
  // Calculate Psirmp and Zrmp
  // .........................
  for (int i = 0; i < NDIAG; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jp = 0; jp < J; jp++)
	    {
	      complex<double> sump1 = - YYY(j,   jp, i);
	      complex<double> sumz1 = - YYY(J+j, jp, i);
	      
	      for (int k = 0; k < nres; k++)
		{
		  sump1 += Psiu(j, k, i) * Pia(k, jp);
		  sumz1 += Zu  (j, k, i) * Pia(k, jp);
		}

	      sump += sump1 * Upsilon[jp];
	      sumz += sumz1 * Upsilon[jp];
	    }

	  Psirmp(j, i) = sump;
	  Zrmp  (j, i) = sumz;
	}
    }

  // ..............................................
  // Calculate Psi_x and Psi_rmp on plasma boundary
  // ..............................................
  for (int i = 0; i <= Nw; i++)
    {
      double theta = tbound[i];
      
      complex<double> sump = complex<double> (0., 0.);
      complex<double> sumx = complex<double> (0., 0.);
      
      for (int j = 0; j < J; j++)
	{
	  sump += Psirmp(j, NDIAG-1) * complex<double> (cos (mpol[j] * theta), sin (mpol[j] * theta));
	  sumx += Psix[j]            * complex<double> (cos (mpol[j] * theta), sin (mpol[j] * theta));
	}
      
      Psirmps[i] = sump;
      Psixs  [i] = sumx;
    }
 }
