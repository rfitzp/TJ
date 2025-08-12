// Matrix.cpp

#include "Vertical.h"

// ##############################################
// Function to return values of coupling matrices
// ##############################################
void Vertical::GetMatrices (double r, 
			    Array<complex<double>,2> AAmmp, Array<complex<double>,2> BBmmp,
			    Array<complex<double>,2> CCmmp, Array<complex<double>,2> DDmmp)
{
  complex<double> Ammp, Bmmp, Cmmp, Dmmp;
  
  for (int j = 0; j < J; j++)
    {
      int m = MPOL[j];

      for (int jp = 0; jp < J; jp++)
	{
	  int mp = MPOL[jp];

	  GetMatrices (r, m, mp, Ammp, Bmmp, Cmmp, Dmmp);

	  AAmmp(j, jp) = Ammp;
	  BBmmp(j, jp) = Bmmp;
	  CCmmp(j, jp) = Cmmp;
	  DDmmp(j, jp) = Dmmp;
	}
    }
}
  
// #####################################################
// Function to return values of coupling matrix elements
// #####################################################
void Vertical::GetMatrices (double r, int m, int mp, 
			    complex<double>& Ammp, complex<double>& Bmmp,
			    complex<double>& Cmmp, complex<double>& Dmmp)
{
  double pp   = Getpp  (r);
  double ppp  = Getppp (r);
  double q    = Getq   (r);
  double s    = Gets   (r);
  double S1   = GetS1  (r);
  double S5   = GetS5  (r);
  double Sigp = GetSigp (r);
  double P1   = GetP1  (r);
  double P2   = GetP2  (r);

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
  double eps2 = epsa*epsa;
  double r2   = r*r;
  double q2   = q*q;
  double mm2  = mm*mm;

  if (m == mp)
    {
      double dmm0 = mm2 + q * P2;
      double dmm2 = mm2 * S5 + (- r2*P1*P1 + q*r*Sigp - r*pp - r2*ppp - 2. * (1. - s) * r*pp + 2.*r*pp*q2 * (-2. + 3.*pp*q2/r)
			       + 2.*Hnp[1]*q2 * (pp + r * ppp - 4. * (1. - s) * pp));

      Ammp = complex<double> (1. + eps2 * (- 0.75*r2 + Hn[1] + S1),       0.);
      Bmmp = complex<double> (0.,                                         0.);
      Cmmp = complex<double> (0.,                                         0.);
      Dmmp = complex<double> (dmm0 + eps2 * dmm2,                         0.);
    }
  else if (mp == m+1)
    {
      Ammp = complex<double> (- epsa * Hnp[1],                                                         0.);
      Bmmp = complex<double> (+ epsa * mmp * (r - pp*q2 + (1. - s) * Hnp[1]),                          0.);
      Cmmp = complex<double> (+ epsa * mm  * (r - pp*q2 + (1. - s) * Hnp[1]),                          0.);
      Dmmp = complex<double> (+ epsa * (pp + r*ppp - (2. - s) * pp) * q2 + epsa*mm*mmp * (r - Hnp[1]), 0.);
    }
  else if (mp == m-1)
    {
      Ammp = complex<double> (- epsa * Hnp[1],                                                         0.);
      Bmmp = complex<double> (- epsa * mmp * (r - pp*q2 + (1. - s) * Hnp[1]),                          0.);
      Cmmp = complex<double> (- epsa * mm  * (r - pp*q2 + (1. - s) * Hnp[1]),                          0.);
      Dmmp = complex<double> (+ epsa * (pp + r*ppp - (2. - s) * pp) * q2 + epsa*mm*mmp * (r - Hnp[1]), 0.);
    }
  else if (mp > m+1)
    {
      int    n  = mp - m;
      double nn = double (n);
      if (n > Ns)
	{
	  Ammp = complex<double> (0., 0.);
	  Bmmp = complex<double> (0., 0.);
	  Cmmp = complex<double> (0., 0.);
	  Dmmp = complex<double> (0., 0.);
	}
      else
	{
	  Ammp = complex<double> (- epsa * Hnp[n],
				  - epsa * Vnp[n]);
	  Bmmp = complex<double> (+ epsa * (mmp /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mmp /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Cmmp = complex<double> (+ epsa * (mm  /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mm  /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Dmmp = complex<double> (- epsa * mm * mmp * Hnp[n],
				  - epsa * mm * mmp * Vnp[n]);
	}
    }
  else if (mp < m-1)
    {
      int    n  = m - mp;
      double nn = double (n);
      if (n > Ns)
	{
	  Ammp = complex<double> (0., 0.);
	  Bmmp = complex<double> (0., 0.);
	  Cmmp = complex<double> (0., 0.);
	  Dmmp = complex<double> (0., 0.);
	}
      else
	{
	  Ammp = complex<double> (- epsa * Hnp[n],
				  + epsa * Vnp[n]);
	  Bmmp = complex<double> (- epsa * (mmp /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mmp /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Cmmp = complex<double> (- epsa * (mm  /nn) * ((1. - s) * Hnp[n] - (nn*nn - 1.) * Hn[n] /r),
				  + epsa * (mm  /nn) * ((1. - s) * Vnp[n] - (nn*nn - 1.) * Vn[n] /r));
	  Dmmp = complex<double> (- epsa * mm * mmp * Hnp[n],
				  + epsa * mm * mmp * Vnp[n]);
	}
    }

  delete[] Hn; delete[] Hnp; delete[] Vn; delete[] Vnp;
 }

