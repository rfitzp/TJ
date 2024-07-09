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
  Pdag.resize (J, J);
  Pinv.resize (J, J);
  Rvac.resize (J, J);
  Rdag.resize (J, J);
  Avac.resize (J, J);
  Hmat.resize (J, J);
  Hdag.resize (J, J);
  Hinv.resize (J, J);

  // .........................
  // Calculate vacuum matrices
  // .........................
  int              neqns = 2*J*J;
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
	Rvac(j, jp) = Y[index]; index++;
      }

  // ..................................
  // Calculate vacuum residual matrix
  // ..................................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> suma = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  {
	    suma += conj (Pvac (jpp, j)) * Rvac (jpp, jp) - conj (Rvac (jpp, j)) * Pvac (jpp, jp);
	  }

	Avac (j, jp) = suma;
      }

  double Aerr = 0.;

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double aval = abs (Avac (j, jp));

	if (aval > Aerr)
	  Aerr = aval;	
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

  SolveLinearSystem (Rdag, Hmat, Pdag);
   
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hdag (j, jp) = conj (Hmat(jp ,j));
   
  double Herr = 0.;

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double hval = abs (Hmat (j, jp) - Hdag (j, jp));

	if (hval > Herr)
	  Herr = hval;	
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

  printf ("Aerr = %11.4e Herr = %10.4e\n", Aerr, Herr);
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
	double Ptorz = ToroidaldPdz (NTOR, MMP, z);

	double Pfac;
	if (MP == 0)
	  {
	    Pfac = sqrt (M_PI) * gsl_sf_gamma (0.5 - ntor) /sqrt (2.);
	  }
	else
	  {
	    Pfac =
	      cos (mmp*M_PI) * sqrt (M_PI) * gsl_sf_gamma (mmp + 0.5 - ntor) * pow (epsa, mmp)
	      /pow (2., mmp - 0.5) /gsl_sf_fact (MMP);
	  }

	complex<double> eik = complex<double> (cos (m * t + mp * eta), - sin (m * t + mp * eta));

	complex<double> Prhs = Pfac * fac * Ptor * eik /2./M_PI;
	complex<double> Rrhs = Pfac * ( (Ptor /2./fac + fac * Ptorz) * R2rz
				       + (set /2./fac - complex<double> (0., 1.) * mp * fac) * Ptor * R2re) * eik /2./M_PI;

	dYdt[index] = Prhs; index++;
	dYdt[index] = Rrhs; index++;
      }
}
