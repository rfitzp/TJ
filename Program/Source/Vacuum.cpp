// Vacuum.cpp

#include "TJ.h"

// #####################################
// Function to calculate vacuum matrices
// #####################################
void TJ::GetVacuum ()
{
  // ...............
  // Allocate memory
  // ...............
  Pvac.resize (J, J);
  Qvac.resize (J, J);
  Rvac.resize (J, J);
  Svac.resize (J, J);
  Avac.resize (J, J);
  Bvac.resize (J, J);
  Cvac.resize (J, J);
  Rdag.resize (J, J);
  Pdag.resize (J, J);
  Hmat.resize (J, J);
  Hdag.resize (J, J);
  Hsym.resize (J, J);
  Imat.resize (J, J);
  Gmat.resize (J, J);
  double* Hn  = new double[Ns+1];
  double* Vn  = new double[Ns+1];
  double* Hnp = new double[Ns+1];
  double* Vnp = new double[Ns+1];

  // .........................
  // Calculate vacuum matrices
  // .........................
  for (int n = 1; n <= Ns; n++)
    {
      Hn [n] = GetHn  (n, 1.);
      Hnp[n] = GetHnp (n, 1.);
    }
  for (int n = 2; n <= Ns; n++)
    {
      Vn [n] = GetVn  (n, 1.);
      Vnp[n] = GetVnp (n, 1.);
    }
  
  double nt2  = ntor*ntor;
  double eps2 = epsa*epsa;
  int    nt   = int (ntor);

  double sum = 0.;
  for (int i = 1; i <= nt; i++)
    sum += 2. /double (2*i - 1);

  double zetan = exp (sum);
  double lnn   = log (8./zetan/epsa);

  for (int jp = 0; jp < J; jp++)
    for (int j = 0; j < J; j++)
      {
	Pvac (j, jp) = complex<double> (0., 0.);
	Qvac (j, jp) = complex<double> (0., 0.);
	Rvac (j, jp) = complex<double> (0., 0.);
	Svac (j, jp) = complex<double> (0., 0.);
      }

  for (int jp = 0; jp < J; jp++)
    {
      int    mp = MPOL[jp];
      double mm = fabs (mpol[jp]);
      double m2 = mm*mm;
	 
      double G0p;
      if (abs (mp) == 1)
	G0p = - 0.25*nt2 - (0.5*nt2 - 0.125) * lnn - 0.5*Hn[1] - 0.25*Hnp[1];
      else
	G0p =      - (nt2 + (mm - 2.) * (mm - 0.75)) /4. /(mm - 1.) - 0.5*Hn[1] - 0.25*Hnp[1];
      double G0m = + (nt2 + (mm + 2.) * (mm + 0.75)) /4. /(mm + 1.) - 0.5*Hn[1] - 0.25*Hnp[1];
      double G3p;
      if (abs (mp) == 1)
	G3p = - 0.25*nt2 + 0.125 + (0.5*nt2 - 0.125) * lnn;
      else
	G3p =      - (  (mm - 2.) * nt2 + 0.25*m2) /4. /(mm - 1.) /mm;
      double G3m = - (- (mm + 2.) * nt2 + 0.25*m2) /4. /(mm + 1.) /mm;
	 
      if (mp == 0)
	{
	  Pvac (jp, jp) = complex<double> (lnn
					   + eps2 * (0.5*nt2 - 0.3125 + (0.25*nt2 + 0.375) * lnn
						     - (0.5*Hn[1] + 0.25*Hnp[1]) * lnn - G1), 0.);
	  Qvac (jp, jp) = complex<double> (+ 1. + eps2 * G0m,                   0.);
	  Rvac (jp, jp) = complex<double> (- 1. + eps2 * 0.5*nt2 * (lnn + 0.5), 0.);
	  Svac (jp, jp) = complex<double> (       eps2 * 0.5*nt2,               0.);

	  if (jp + 1 < J)
	    {
	      Pvac (jp+1, jp) = complex<double> (epsa * (0.25*lnn - 0.25 + 0.5*Hn[1]),              0.);
	      Qvac (jp+1, jp) = complex<double> (epsa * (0.25),                                     0.);
	      Rvac (jp+1, jp) = complex<double> (epsa * (0.25*lnn + 0.5  - 0.5*Hn[1] - 0.5*Hnp[1]), 0.);
	      Svac (jp+1, jp) = complex<double> (epsa * (0.25),                                     0.);
	    }

	  if (jp - 1 > -1)
	    {
	      Pvac (jp-1, jp) = complex<double> (epsa * (0.25*lnn - 0.25 + 0.5*Hn[1]),              0.);
	      Qvac (jp-1, jp) = complex<double> (epsa * (0.25),                                     0.);
	      Rvac (jp-1, jp) = complex<double> (epsa * (0.25*lnn + 0.5  - 0.5*Hn[1] - 0.5*Hnp[1]), 0.);
	      Svac (jp-1, jp) = complex<double> (epsa * (0.25),                                     0.);
	    }

	  for (int k = 2; k <= Ns; k++)
	    {
	      if (jp + k < J)
		{
		  Pvac (jp+k, jp) = complex<double> (+ epsa * 0.5*Hn[k],              - epsa * 0.5*Vn[k]);
		  Qvac (jp+k, jp) = complex<double> (  0.,                              0.);
		  Rvac (jp+k, jp) = complex<double> (- epsa * 0.5 * (Hnp[k] + Hn[k]), + epsa * 0.5 * (Vnp[k] + Vn[k]));
		  Svac (jp+k, jp) = complex<double> (  0.,                              0.);
		}

	       if (jp - k > -1)
		 {
		   Pvac (jp-k, jp) = complex<double> (+ epsa * 0.5*Hn[k],              + epsa * 0.5*Vn[k]);
		   Qvac (jp-k, jp) = complex<double> (  0.,                              0.);
		   Rvac (jp-k, jp) = complex<double> (- epsa * 0.5 * (Hnp[k] + Hn[k]), - epsa * 0.5 * (Vnp[k] + Vn[k]));
		   Svac (jp-k, jp) = complex<double> (  0.,                              0.);
		}
  	    }
	}
      else
	{
	  int    sig;
	  double sigma;
	  if (mp > 0)
	    {
	      sig   = 1;
	      sigma = 1.;
	    }
	  else
	    {
	      sig   = - 1;
	      sigma = - 1.;
	    }

	  Pvac (jp, jp) = complex<double> (      + 1. + eps2 * (G0p - mm *  G1              - m2 * G2),  0.);
	  Qvac (jp, jp) = complex<double> (      + 1. + eps2 * (G0m + mm *  G1              - m2 * G2),  0.);
	  Rvac (jp, jp) = complex<double> (mm * (- 1. - eps2 * (G3p - mm * (G1 + 0.5*Hn[1]) - m2 * G2)), 0.);
	  Svac (jp, jp) = complex<double> (mm * (+ 1. + eps2 * (G3m + mm * (G1 + 0.5*Hn[1]) - m2 * G2)), 0.);

	  if (jp + sig < J && jp + sig > -1)
	    {
	      Pvac (jp+sig, jp) = complex<double> (epsa             * (0.25 - 0.5 *mm + mm * (Hn[1] + 0.5*Hnp[1])),   0.);
	      Qvac (jp+sig, jp) = complex<double> (epsa             * (0.25 + 0.5 *mm * Hnp[1]),                      0.);
	      Rvac (jp+sig, jp) = complex<double> (epsa * (1. + mm) * (0.25 + 0.5 *mm - mm * (Hn[1] + 0.5*Hnp[1])),   0.);
	      Svac (jp+sig, jp) = complex<double> (epsa             * (0.25 - 0.25*mm + 0.5*mm * (1. + mm) * Hnp[1]), 0.);
	    }

	  if (jp - sig < J && jp - sig > -1)
	    {
	      Pvac (jp-sig, jp) = complex<double> (epsa             * (0.25 - 0.5 *mm * Hnp[1]),                      0.);
	      Qvac (jp-sig, jp) = complex<double> (epsa             * (0.25 + 0.5 *mm - mm * (Hn[1] + 0.5*Hnp[1])),   0.);
	      Rvac (jp-sig, jp) = complex<double> (epsa             * (0.25 + 0.25*mm - 0.5*mm * (1. - mm) * Hnp[1]), 0.);
	      Svac (jp-sig, jp) = complex<double> (epsa * (1. - mm) * (0.25 - 0.5 *mm + mm * (Hn[1] + 0.5*Hnp[1])),   0.);
	    }

	  for (int k = 2; k <= Ns; k++)
	    {
	      double kk = double (k);
	      
	      if (jp + sig*k < J && jp + sig*k > -1)
		{
		  Pvac (jp+sig*k, jp) = complex<double> (epsa * (mm/2./kk)                * (+ Hnp[k] + (kk + 1.) * Hn[k]),
							 epsa * (mm/2./kk)                * (- Vnp[k] - (kk + 1.) * Vn[k]) * sigma);
		  Qvac (jp+sig*k, jp) = complex<double> (epsa * (mm/2./kk)                * (+ Hnp[k] - (kk - 1.) * Hn[k]),
							 epsa * (mm/2./kk)                * (- Vnp[k] + (kk - 1.) * Vn[k]) * sigma);
		  Rvac (jp+sig*k, jp) = complex<double> (epsa * (mm/2.   ) * (1. + mm/kk) * (- Hnp[k] - (kk + 1.) * Hn[k]),
							 epsa * (mm/2.   ) * (1. + mm/kk) * (+ Vnp[k] + (kk + 1.) * Vn[k]) * sigma);
		  Svac (jp+sig*k, jp) = complex<double> (epsa * (mm/2.   ) * (1. + mm/kk) * (+ Hnp[k] - (kk - 1.) * Hn[k]),
							 epsa * (mm/2.   ) * (1. + mm/kk) * (- Vnp[k] + (kk - 1.) * Vn[k]) * sigma);
		}
	      
	      if (jp - sig*k < J && jp - sig*k > -1)
		{
		  Pvac (jp-sig*k, jp) = complex<double> (epsa * (mm/2./kk)                * (- Hnp[k] + (kk - 1.) * Hn[k]),
							 epsa * (mm/2./kk)                * (- Vnp[k] + (kk - 1.) * Vn[k]) * sigma);
		  Qvac (jp-sig*k, jp) = complex<double> (epsa * (mm/2./kk)                * (- Hnp[k] - (kk + 1.) * Hn[k]),
							 epsa * (mm/2./kk)                * (- Vnp[k] - (kk + 1.) * Vn[k]) * sigma);
		  Rvac (jp-sig*k, jp) = complex<double> (epsa * (mm/2.   ) * (1. - mm/kk) * (- Hnp[k] + (kk - 1.) * Hn[k]),
							 epsa * (mm/2.   ) * (1. - mm/kk) * (- Vnp[k] + (kk - 1.) * Vn[k]) * sigma);
		  Svac (jp-sig*k, jp) = complex<double> (epsa * (mm/2.   ) * (1. - mm/kk) * (+ Hnp[k] + (kk + 1.) * Hn[k]),
							 epsa * (mm/2.   ) * (1. - mm/kk) * (+ Vnp[k] + (kk + 1.) * Vn[k]) * sigma);
		}
	    }
	}
    }

  delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;

  // ...........................
  // Renormalize vacuum matrices
  // ...........................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	int    mp = MPOL[jp];
	double mm = fabs (mpol[jp]);

	if (mp != 0)
	  {
	    Pvac (j, jp) /= mm;
	    Qvac (j, jp) /= mm;
	    Rvac (j, jp) /= mm;
	    Svac (j, jp) /= mm;
	  }
      }

  // ..................................
  // Calculate vacuum residual matrices
  // ..................................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> suma = complex<double> (0., 0.);
	complex<double> sumb = complex<double> (0., 0.);
	complex<double> sumc = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  {
	    suma += conj (Pvac (jpp, j)) * Rvac (jpp, jp) - conj (Rvac (jpp, j)) * Pvac (jpp, jp);
	    sumb += conj (Qvac (jpp, j)) * Svac (jpp, jp) - conj (Svac (jpp, j)) * Qvac (jpp, jp);
	    sumc += conj (Pvac (jpp, j)) * Svac (jpp, jp) - conj (Rvac (jpp, j)) * Qvac (jpp, jp);
	  }

	Avac (j, jp) = suma;
	Bvac (j, jp) = sumb;
	Cvac (j, jp) = sumc;

	if (j == jp)
	  {
	    if (MPOL[j] == 0)
	      Cvac (j, jp) -= complex<double> (1., 0.);
	    else
	      Cvac (j, jp) -= complex<double> (2./fabs (mpol[j]), 0.);
	  }
      }

  // ..................
  // Calculate H-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Pdag (j, jp) = conj (Pvac (jp, j));
	Rdag (j, jp) = conj (Rvac (jp, j));
      }
  
  SolveLinearSystem (Pdag, Hmat, Rdag);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hdag (j, jp) = conj (Hmat(jp ,j));
   
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hsym(j, jp) = 0.5 * (Hmat(j, jp) + Hdag(j, jp));

  // ..................
  // Calculate G-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = Svac (j, jp);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Hmat (j, jpp) * Qvac (jpp, jp);

	Gmat (j, jp) = sum;
      }
}
 
