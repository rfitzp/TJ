// Matrix.cpp

#include "TJ.h"

// ##############################################
// Function to return values of coupling matrices
// ##############################################
void TJ::GetMatrices (double r,
		      Array<complex<double>,2> LLmmp, Array<complex<double>,2> MMmmp,
		      Array<complex<double>,2> NNmmp, Array<complex<double>,2> PPmmp)
{
  complex<double> Lmmp, Mmmp, Nmmp, Pmmp;
  
  for (int j = 0; j < J; j++)
    {
      int m = MPOL[j];

      for (int jp = 0; jp < J; jp++)
	{
	  int mp = MPOL[jp];

	  GetMatrices (r, m, mp, Lmmp, Mmmp, Nmmp, Pmmp);

	  LLmmp(j, jp) = Lmmp;
	  MMmmp(j, jp) = Mmmp;
	  NNmmp(j, jp) = Nmmp;
	  PPmmp(j, jp) = Pmmp;
	}
    }
}
  
// #####################################################
// Function to return values of coupling matrix elements
// #####################################################
void TJ::GetMatrices (double r, int m, int mp,
		      complex<double>& Lmmp, complex<double>& Mmmp,
		      complex<double>& Nmmp, complex<double>& Pmmp)
{
  double pp  = Getpp  (r);
  double ppp = Getppp (r);
  double q   = Getq   (r);
  double s   = Gets   (r);
  double s2  = Gets2  (r);
  double S1  = GetS1  (r);
  double P1  = GetP1  (r);
  double P2  = GetP2  (r);
  double P3  = GetP3  (r);

  double* Hn  = new double[Ns+1];
  double* Hnp = new double[Ns+1];
  double* Vn  = new double[Ns+1];
  double* Vnp = new double[Ns+1];

  for (int n = 1; n <= Ns; n++)
    {
      Hn[n]  = GetHn  (n, r);
      Hnp[n] = GetHnp (n, r);
    }
  for (int n = 2; n <= Ns; n++)
    {
      Vn[n]  = GetVn  (n, r);
      Vnp[n] = GetVnp (n, r);
    }

  double mm   = double (m);
  double mmp  = double (mp);
  double mnq  = mm  - ntor*q;
  double mpnq = mmp - ntor*q;
  double eps2 = epsa*epsa;
  double r2   = r*r;
  double q2   = q*q;
  double q3   = q2*q;
  double mm2  = mm*mm;
  double nt2  = ntor*ntor;
  
  if (m == 0 && mp == 0)
    {
      Lmmp = complex<double> (eps2*nt2*r2,                          0.);
      Mmmp = complex<double> (0.,                                   0.);
      Nmmp = complex<double> (0.,                                   0.);
      Pmmp = complex<double> (q2 * (nt2 - 2.*P1 - P2 - ppp + pp/r), 0.);
    }
  else if (m == mp)
    {
      double pmm0 = mnq*mnq + mnq * (q /mm) * P2;
      double pmm2 = mnq*mnq * (1.75*r2 - Hn[1] - 3.*r*Hnp[1] + S1
			       + (1./mm2) * ((ntor /mm) * r2 * (2.*P1 + P2) - r2*P1*P1 - r*pp - r2*ppp))
	- (mnq /mm) * P3 + 2*r*pp*(1. - q2);
      
      Lmmp = complex<double> (mm2 + eps2 * mm2 * (- 0.75*r2 + Hn[1] + S1) + eps2*nt2*r2, 0.);
      Mmmp = complex<double> (0.,                                                        0.);
      Nmmp = complex<double> (0.,                                                        0.);
      Pmmp = complex<double> (pmm0 + eps2 * pmm2,                                        0.);
    }
  else if (mp == m+1)
    {
      Lmmp = complex<double> (- epsa * mm*mmp   * Hnp[1],                                            0.);
      Mmmp = complex<double> (- epsa * mm*mnq   * pp*q2 + epsa * mm*mpnq  * (r + (1. - s) * Hnp[1]), 0.);
      Nmmp = complex<double> (- epsa * mmp*mpnq * pp*q2 + epsa * mmp*mnq  * (r + (1. - s) * Hnp[1]), 0.);
      Pmmp = complex<double> (- epsa * (1. + s) * pp*q2 + epsa * mnq*mpnq * (r - Hnp[1]),            0.);

      if (mp == 0)
	{
	  Mmmp += complex<double> (  epsa * (2. - s) * Hnp[1],                                                   0.);
	  Pmmp += complex<double> (- epsa * (2. - s) * (- ntor*q3*pp + (1. + ntor*q) * (r + (1. - s) * Hnp[1])), 0.);
	}
      if (m == 0)
	{
	  Nmmp += complex<double> (  epsa * (2. - s) * Hnp[1],                                                   0.);
	  Pmmp += complex<double> (- epsa * (2. - s) * (+ ntor*q3*pp + (1. - ntor*q) * (r + (1. - s) * Hnp[1])), 0.);
	}
    }
  else if (mp == m-1)
    {
      Lmmp = complex<double> (- epsa * mm*mmp   * Hnp[1],                                            0.);
      Mmmp = complex<double> (+ epsa * mm*mnq   * pp*q2 - epsa * mm*mpnq  * (r + (1. - s) * Hnp[1]), 0.);
      Nmmp = complex<double> (+ epsa * mmp*mpnq * pp*q2 - epsa * mmp*mnq  * (r + (1. - s) * Hnp[1]), 0.);
      Pmmp = complex<double> (- epsa * (1. + s) * pp*q2 + epsa * mnq*mpnq * (r - Hnp[1]),            0.);

      if (mp == 0)
	{
	  Mmmp += complex<double> (- epsa * (2. - s) * Hnp[1],                                                   0.);
	  Pmmp += complex<double> (- epsa * (2. - s) * (+ ntor*q3*pp + (1. - ntor*q) * (r + (1. - s) * Hnp[1])), 0.);
	}
      if (m == 0)
	{
	  Nmmp += complex<double> (- epsa * (2. - s) * Hnp[1],                                                   0.);
	  Pmmp += complex<double> (- epsa * (2. - s) * (- ntor*q3*pp + (1. + ntor*q) * (r + (1. - s) * Hnp[1])), 0.);
	}
    }
  else if (mp > m+1)
    {
      int    n  = mp - m;
      double nn = double (n);
      if (n > Ns)
	{
	  Lmmp = complex<double> (0., 0.);
	  Mmmp = complex<double> (0., 0.);
	  Nmmp = complex<double> (0., 0.);
	  Pmmp = complex<double> (0., 0.);
	}
      else
	{
	  Lmmp = complex<double> (- epsa * mm*mmp * Hnp[n],
				  - epsa * mm*mmp * Vnp[n]);
	  Mmmp = complex<double> (+ epsa * (mm /nn)  * mpnq * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mm /nn)  * mpnq * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Nmmp = complex<double> (+ epsa * (mmp /nn) * mnq  * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mmp /nn) * mnq  * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Pmmp = complex<double> (- epsa * mnq*mpnq * Hnp[n],
				  - epsa * mnq*mpnq * Vnp[n]);

	  if (mp == 0)
	    {
	      Mmmp += complex<double> (+ epsa * (2. - s) * nn * Hnp[n],
				       + epsa * (2. - s) * nn * Vnp[n]);
	      Pmmp += complex<double> (- epsa * (2. - s) * ((nn + ntor*q) /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				       - epsa * (2. - s) * ((nn + ntor*q) /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	    }
	  if (m == 0)
	    {
	      Nmmp += complex<double> (+ epsa * (2. - s) * nn * Hnp[n],
				       + epsa * (2. - s) * nn * Vnp[n]);
	      Pmmp += complex<double> (- epsa * (2. - s) * ((nn - ntor*q) /nn)* ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				       - epsa * (2. - s) * ((nn - ntor*q) /nn)* ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	    }
	}
    }
  else if (mp < m-1)
    {
      int    n  = m - mp;
      double nn = double (n);
      if (n > Ns)
	{
	  Lmmp = complex<double> (0., 0.);
	  Mmmp = complex<double> (0., 0.);
	  Nmmp = complex<double> (0., 0.);
	  Pmmp = complex<double> (0., 0.);
	}
      else
	{
	  Lmmp = complex<double> (- epsa * mm*mmp * Hnp[n],
				  + epsa * mm*mmp * Vnp[n]);
	  Mmmp = complex<double> (- epsa * (mm /nn)  * mpnq * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mm /nn)  * mpnq * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Nmmp = complex<double> (- epsa * (mmp /nn) * mnq  * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mmp /nn) * mnq  * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Pmmp = complex<double> (- epsa * mnq*mpnq * Hnp[n],
				  + epsa * mnq*mpnq * Vnp[n]);

	  if (mp == 0)
	    {
	      Mmmp += complex<double> (- epsa * (2. - s) * nn * Hnp[n],
				       + epsa * (2. - s) * nn * Vnp[n]);
	      Pmmp += complex<double> (- epsa * (2. - s) * ((nn - ntor*q) /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				       + epsa * (2. - s) * ((nn - ntor*q) /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	    }
	  if (m == 0)
	    {
	      Nmmp += complex<double> (- epsa * (2. - s) * nn * Hnp[n],
				       + epsa * (2. - s) * nn * Vnp[n]);
	      Pmmp += complex<double> (- epsa * (2. - s) * ((nn + ntor*q) /nn)* ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				       + epsa * (2. - s) * ((nn + ntor*q) /nn)* ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	    }
	}
    }

  delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;
}
