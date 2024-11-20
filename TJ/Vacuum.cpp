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
  Pdag.resize (J, J);
  Rdag.resize (J, J);
  Qdag.resize (J, J);
  Amat.resize (J, J);
  Aher.resize (J, J);
  Aant.resize (J, J);
  Bmat.resize (J, J);
  Bher.resize (J, J);
  Bant.resize (J, J);
  Imat.resize (J, J);

  Cmat.resize (J, J);
  Cdag.resize (J, J);
  Hmat.resize (J, J);

  PQmat.resize (J, J);
  RSmat.resize (J, J);
  Dmat .resize (J, J);
  Ddag .resize (J, J);
  Gmat .resize (J, J);

  Psix = new complex<double>[J];
  Xi   = new complex<double>[J];

  // .........................
  // Calculate vacuum matrices
  // .........................
  int              neqns = 4*J*J;
  double           h, t_err, t;
  int              rept;
  complex<double>* Y   = new complex<double>[neqns];
  complex<double>* err = new complex<double>[neqns];
  rhs_chooser          = 0;
  
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
	Qvac(j, jp) = Y[index]; index++;
	Svac(j, jp) = Y[index]; index++;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Pdag (j, jp) = conj (Pvac (jp, j));
	Rdag (j, jp) = conj (Rvac (jp, j));
	Qdag (j, jp) = conj (Qvac (jp, j));
      }

  // ..................
  // Calculate A-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Pdag (j, jpp) * Rvac (jpp, jp);

	Amat (j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Aher (j, jp) = 0.5 * (Amat (j, jp) + conj (Amat (jp, j)));
	Aant (j, jp) = 0.5 * (Amat (j, jp) - conj (Amat (jp, j)));
      }
  
  // ............................
  // Calculate A-matrix residuals
  // ............................
  double Ahmax = 0., Aamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ahval = abs (Aher (j, jp));
	double aaval = abs (Aant (j, jp));

	if (ahval > Ahmax)
	  Ahmax = ahval;
	if (aaval > Aamax)
	  Aamax = aaval;	
      }

  // ..................
  // Calculate B-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Qdag (j, jpp) * Svac (jpp, jp);

	Bmat (j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Bher (j, jp) = 0.5 * (Bmat (j, jp) + conj (Bmat (jp, j)));
	Bant (j, jp) = 0.5 * (Bmat (j, jp) - conj (Bmat (jp, j)));
      }
  
  // ............................
  // Calculate B-matrix residuals
  // ............................
  double Bhmax = 0., Bamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double bhval = abs (Bher (j, jp));
	double baval = abs (Bant (j, jp));

	if (bhval > Bhmax)
	  Bhmax = bhval;
	if (baval > Bamax)
	  Bamax = baval;	
      }

  // ..................
  // Calculate I-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Pdag (j, jpp) * Svac (jpp, jp) - Rdag (j, jpp) * Qvac (jpp, jp);

	Imat (j, jp) = sum;
      }
  
  // ............................
  // Calculate I-matrix residuals
  // ............................
  double Imax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ival;
	if (j == jp)
	  ival = abs (Imat (j, jp) - 1.);
	else
	  ival = abs (Imat (j, jp));

	if (ival > Imax)
	  Imax = ival;
      }
  
  printf ("Vacuum Hermitian test residual: %10.4e %10.4e %10.4e\n", Aamax/Ahmax, Bamax/Bhmax, Imax);

  // ..................
  // Calculate H-matrix
  // ..................
  SolveLinearSystemTranspose (Rvac, Cmat, Pvac);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Cdag (j, jp) = conj (Cmat (jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hmat (j, jp) = 0.5 * Cmat (j, jp) + 0.5 * Cdag (j, jp);

  // ..................
  // Calculate G-matrix
  // ..................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	PQmat (j, jp) = Pvac (j, jp);
	RSmat (j, jp) = Rvac (j, jp);

	for (int jpp = 0; jpp < J; jpp++)
	  {
	    PQmat (j, jp) += Qvac (j, jpp) * Iw (jpp, jp);
	    RSmat (j, jp) += Svac (j, jpp) * Iw (jpp, jp);
	  }
      }
  SolveLinearSystemTranspose (RSmat, Dmat, PQmat);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Ddag (j, jp) = conj (Dmat (jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Gmat (j, jp) = 0.5 * Dmat (j, jp) + 0.5 * Ddag (j, jp);
    
  // ......................................................
  // Calculate poloidal harmonics of RMP at plasma boundary
  // ......................................................
  complex<double>* Y1   = new complex<double>[J];
  complex<double>* err1 = new complex<double>[J];
  rhs_chooser           = 1;
  
  for (int i = 0; i < J; i++)
    Y1[i] = complex<double> (0., 0.);
  
  t     = 0.;
  h     = h0;
  count = 0;
  
  do
    {
      CashKarp45Adaptive1 (J, t, Y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (t < 2.*M_PI - h);
  CashKarp45Fixed1 (J, t, Y1, err1, 2.*M_PI - t);

  for (int i = 0; i < J; i++)
    Psix[i] = Y1[i];
  
  // ...................
  // Calculate Xi-vector
  // ...................
  for (int j = 0; j < J; j++)
    {
      complex<double> sum = complex<double> (0., 0.);

      for (int jp = 0; jp < J; jp++)
	sum += Hmat(j, jp) * Psix[jp];

      Xi[j] = sum;
    }

  delete[] Y; delete[] err; delete[] Y1; delete[] err1;
}
 
// ####################################################
// Function to evaluate right-hand sides of vacuum odes
// ####################################################
void TJ::Rhs1 (double t, complex<double>* Y, complex<double>* dYdt)
{
  if (rhs_chooser == 0)
    {
      // ...................................................
      // Right-hand sides for calculation of vacuum matrices
      // ...................................................
      
      int index = 0;
      for (int j = 0; j < J; j++)
	for (int jp = 0; jp < J; jp++)
	  {
	    int    MP  = MPOL[jp];
	    int    MMP = abs (MP);
	    
	    double m   = mpol[j];
	    double mp  = mpol[jp];
	    
	    double R    = gsl_spline_eval (Rbspline,  t, Rbacc);
	    double Z    = gsl_spline_eval (Zbspline,  t, Zbacc);
	    double R2rz = gsl_spline_eval (Rrzspline, t, Rrzacc);
	    double R2re = gsl_spline_eval (Rrespline, t, Rreacc);
	    
	    double z   = GetCoshMu (R, Z);
	    double eta = GetEta    (R, Z);
	    double cet = cos (eta);
	    double set = sin (eta);
	    
	    double fac = sqrt (z - cet);
	    
	    double Ptor  = NormToroidalP    (NTOR, MMP, z);
	    double Ptorz = NormToroidaldPdz (NTOR, MMP, z);
	    double Qtor  = NormToroidalQ    (NTOR, MMP, z);
	    double Qtorz = NormToroidaldQdz (NTOR, MMP, z);
	    
	    complex<double> eik = complex<double> (cos (m * t + mp * eta), - sin (m * t + mp * eta));

	    complex<double> Prhs = fac * Ptor * eik /2./M_PI;
	    complex<double> Rrhs = (  (Ptor /2./fac + fac * Ptorz) * R2rz
                                    + (set /2./fac - complex<double> (0., 1.) * mp * fac) * Ptor * R2re) * eik /2./M_PI;
	    complex<double> Qrhs = fac * Qtor * eik /2./M_PI;
	    complex<double> Srhs = (  (Qtor /2./fac + fac * Qtorz) * R2rz
                                    + (set /2./fac - complex<double> (0., 1.) * mp * fac) * Qtor * R2re) * eik /2./M_PI;
	    
	    dYdt[index] = Prhs; index++;
	    dYdt[index] = Rrhs; index++;
	    dYdt[index] = Qrhs; index++;
	    dYdt[index] = Srhs; index++;
	  }
    }
  else if (rhs_chooser == 1)
    {
      // ................................................................................
      // Right-hand sides for calculation of poloidal harmonics of RMP at plasma boundary
      // ................................................................................

      double R = gsl_spline_eval (Rbspline, t, Rbacc);
      double Z = gsl_spline_eval (Zbspline, t, Zbacc);

      double Aphi;
      for (int i = 0; i < ncoil; i++)
	Aphi += Icoil[i] * GetG (R, Z, Rcoil[i], Zcoil[i]);

      for (int j = 0; j < J; j++)
	dYdt[j] = mpol[j] * R * Aphi * complex<double> (cos (mpol[j] * t), - sin (mpol[j] * t)) /2./M_PI;
    }
}

// ###############################################
// Function to evaluate G for RMP coil calculation
// ###############################################
double TJ::GetG (double R, double Z, double Rp, double Zp)
{
  double eta  = Geteta (R, Z, Rp, Zp);
  double ceta = cosh (eta);

  double fun = (ntor - 0.5) * ToroidalP (NTOR - 1, 0, ceta) + ToroidalP (NTOR + 1, 0, ceta) /(ntor + 0.5);

  fun *= cos ((ntor+1.) * M_PI) * sqrt (M_PI) * R * Rp /4. /gsl_sf_gamma (ntor + 0.5);
  fun *= sqrt (ceta /(R*R + Rp*Rp + (Z - Zp) * (Z - Zp)));

  return fun;
}

// #################################################
// Function to evaluate eta for RMP coil calculation
// #################################################
double TJ::Geteta (double R, double Z, double Rp, double Zp)
{
  double tanhe = 2. * R * Rp /(R*R + Rp*Rp + (Z - Zp) * (Z - Zp));

  return atanh (tanhe);
}


