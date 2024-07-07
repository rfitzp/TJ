// Vacuum.cpp

#include "TJ.h"

// #####################################
// Function to calculate vacuum matrices
// #####################################
void TJ::GetVacuum ()
{
  printf ("Calculating vacuum data:\n");
  
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
  Pinv.resize (J, J);
  Hinv.resize (J, J);
  Hmat.resize (J, J);
  Hdag.resize (J, J);
  Hsym.resize (J, J);
  Gmat.resize (J, J);

  // .........................
  // Calculate vacuum matrices
  // .........................
  if (TVAC)
    {
      int              neqns = 4*J*J;
      double           h, t_err, t;
      int              rept;
      complex<double>* Y   = new complex<double>[neqns];
      complex<double>* err = new complex<double>[neqns];

      for (int i = 0; i < neqns; i++)
	Y[i] = complex<double> (0., 0.);

      t     = 0.;
      h     = h0;
      count = 0;

      do
	{
	  CashKarp45Adaptive1 (neqns, t, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (t < 2.*M_PI - h);
      CashKarp45Fixed1 (neqns, t, Y, err, 2.*M_PI - t);

      int index = 0;
      for (int j = 0; j < J; j++)
	for (int jp = 0; jp < J; jp++)
	  {
	    Pvac(j, jp) = Y[index]; index++;
	    Qvac(j, jp) = Y[index]; index++;
	    Rvac(j, jp) = Y[index]; index++;
	    Svac(j, jp) = Y[index]; index++;
	  }
     }
  else
    {
      double* Hn  = new double[Ns+1];
      double* Vn  = new double[Ns+1];
      double* Hnp = new double[Ns+1];
      double* Vnp = new double[Ns+1];

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

  double Aerr = 0.;
  double Berr = 0.;
  double Cerr = 0.;

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double aval = abs (Avac (j, jp));
	double bval = abs (Bvac (j, jp));
	double cval = abs (Cvac (j, jp));

	if (aval > Aerr)
	  Aerr = aval;	
	if (bval > Berr)
	  Berr = bval;
	if (cval > Cerr)
	  Cerr = cval;
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
 
  if (SYMM)
    {
      InvertMatrix (Pvac, Pinv);

      for (int j = 0; j < J; j++)
	for (int jp = 0; jp < J; jp++)
	  {
	    complex<double> sum = complex<double> (0., 0.);

	    for (int jpp = 0; jpp < J; jpp++)
	      Hinv (j, jp) += 0.5 * (Rvac (j, jpp) * Pinv (jpp, jp) + conj (Rvac (jp, jpp) * Pinv (jpp, j)));
	  }

      InvertMatrix (Hinv, Hmat);
    }
  else
    {
      SolveLinearSystem (Rdag, Hmat, Pdag);
    }
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hdag (j, jp) = conj (Hmat(jp ,j));
   
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hsym(j, jp) = 0.5 * (Hmat(j, jp) + Hdag(j, jp));

  double Herr = 0.;

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double hval = abs (Hmat (j, jp) - Hdag (j, jp));

	if (hval > Herr)
	  Herr = hval;	
      }

  printf ("Aerr = %11.4e Berr = %10.4e Cerr = %10.4e Herr = %10.4e\n",
	  Aerr, Berr, Cerr, Herr);

  // ..................
  // Calculate G-matrix
  // ..................
  Array<complex<double>,2> Imat(J, J);
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	if (j == jp)
	  {
	    if (MPOL[j] == 0)
	      Imat (j, jp) = complex<double> (1., 0.);
	    else
	      Imat (j, jp) = complex<double> (2./fabs (mpol[j]), 0.);
	  }
	else
	  Imat (j, jp) = complex<double> (0., 0.);
      }
  
  SolveLinearSystem (Rdag, Gmat, Imat);
}
 
// ####################################################
// Function to evaluate right-hand sides of vacuum odes
// ####################################################
void TJ::Rhs1 (double t, complex<double>* Y, complex<double>* dYdt)
{
  int index = 0;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	int    M   = MPOL[j];
	int    MP  = MPOL[jp];
	int    MM  = abs (M);
	int    MMP = abs (MP);
       
	double m   = mpol[j];
	double mp  = mpol[jp];
	double mm  = fabs (m);
	double mmp = fabs (mp);
	
	double R    = gsl_spline_eval (Rbspline,  t, Rbacc);
	double Z    = gsl_spline_eval (Zbspline,  t, Zbacc);
	double R2rz = gsl_spline_eval (Rrzspline, t, Rrzacc);
	double R2re = gsl_spline_eval (Rrespline, t, Rreacc);

	double z   = GetCoshMu (R, Z);
	double eta = GetEta    (R, Z);
	double cet = cos (eta);
	double set = sin (eta);

	double fac = sqrt (z - cet);

	double Ptor  = ToroidalP    (NTOR, MMP, z);
	double Qtor  = ToroidalQ    (NTOR, MMP, z);
	double Ptorz = ToroidaldPdz (NTOR, MMP, z);
	double Qtorz = ToroidaldQdz (NTOR, MMP, z);

	double Pfac, Qfac;
	if (MP == 0)
	  {
	    Pfac = sqrt (M_PI) * gsl_sf_gamma (0.5 - ntor) /sqrt (2.);
	    Qfac = sqrt (2.) /sqrt (M_PI) /gsl_sf_gamma (0.5 + ntor);
	  }
	else
	  {
	    Pfac =
	      cos (mmp*M_PI) * sqrt (M_PI) * gsl_sf_gamma (mmp + 0.5 - ntor) * pow (epsa, mmp)
	      /pow (2., mmp - 0.5) /gsl_sf_fact (MMP);
	    Qfac =
	      cos (mmp*M_PI) * pow (2., mmp + 0.5) * gsl_sf_fact (MMP - 1)
	      /sqrt (M_PI) /gsl_sf_gamma (mmp + 0.5 + ntor) /pow (epsa, mmp);
	  }

	complex<double> eik = complex<double> (cos (m * t + mp * eta), - sin (m * t + mp * eta));

	complex<double> Prhs = Pfac * fac * Ptor * eik /2./M_PI;
	complex<double> Qrhs = Qfac * fac * Qtor * eik /2./M_PI; 
	complex<double> Rrhs = Pfac * ( (Ptor /2./fac + fac * Ptorz) * R2rz
				       + (set /2./fac - complex<double> (0., 1.) * mp * fac) * Ptor * R2re) * eik /2./M_PI;
	complex<double> Srhs = Qfac * ( (Qtor /2./fac + fac * Qtorz) * R2rz
				       + (set /2./fac - complex<double> (0., 1.) * mp * fac) * Qtor * R2re) * eik /2./M_PI;

	dYdt[index] = Prhs; index++;
	dYdt[index] = Qrhs; index++;
	dYdt[index] = Rrhs; index++;
	dYdt[index] = Srhs; index++;
      }
}
