// Vacuum.cpp

#include "TJ.h"

// ##############################################
// Function to calculate vacuum boundary matrices
// ##############################################
void TJ::GetVacuumBoundary ()
{
  printf ("Calculating vacuum boundary data:\n");
  
  // ...............
  // Allocate memory
  // ...............
  Pvac.resize(J, J);
  Qvac.resize(J, J);
  Rvac.resize(J, J);
  Svac.resize(J, J);

  Pdag.resize(J, J);
  Qdag.resize(J, J);
  Rdag.resize(J, J);
  Sdag.resize(J, J);

  PRmat.resize(J, J);
  PRher.resize(J, J);
  PRant.resize(J, J);

  QSmat.resize(J, J);
  QSher.resize(J, J);
  QSant.resize(J, J);

  PSmat.resize(J, J);

  QPmat.resize(J, J);
  QPher.resize(J, J);
  QPant.resize(J, J);

  RSmat.resize(J, J);
  RSher.resize(J, J);
  RSant.resize(J, J);

  SPmat.resize(J, J);

  RPmat.resize(J, J);
  RPdag.resize(J, J);
  Hmat.resize (J, J);

  iRPmat.resize(J, J);
  iRPdag.resize(J, J);
  iHmat .resize(J, J);

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
	Pdag(j, jp) = conj (Pvac(jp, j));
	Qdag(j, jp) = conj (Qvac(jp, j));
	Rdag(j, jp) = conj (Rvac(jp, j));
	Sdag(j, jp) = conj (Svac(jp, j));
      }
  
  // ...................
  // Calculate PR-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Pdag(j, jpp) * Rvac(jpp, jp);

	PRmat(j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	PRher(j, jp) = 0.5 * (PRmat(j, jp) + conj (PRmat(jp, j)));
	PRant(j, jp) = 0.5 * (PRmat(j, jp) - conj (PRmat(jp, j)));
      }
  
  // ............................
  // Calculate PR-matrix residual
  // ............................
  double Ahmax = 0., Aamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ahval = abs (PRher(j, jp));
	double aaval = abs (PRant(j, jp));

	if (ahval > Ahmax)
	  Ahmax = ahval;
	if (aaval > Aamax)
	  Aamax = aaval;	
      }

  // ...................
  // Calculate QS-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Qdag (j, jpp) * Svac (jpp, jp);

	QSmat (j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	QSher(j, jp) = 0.5 * (QSmat(j, jp) + conj (QSmat(jp, j)));
	QSant(j, jp) = 0.5 * (QSmat(j, jp) - conj (QSmat(jp, j)));
      }
  
  // ............................
  // Calculate QS-matrix residual
  // ............................
  double Bhmax = 0., Bamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double bhval = abs (QSher(j, jp));
	double baval = abs (QSant(j, jp));

	if (bhval > Bhmax)
	  Bhmax = bhval;
	if (baval > Bamax)
	  Bamax = baval;	
      }

  // ...................
  // Calculate PS-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Pdag(j, jpp) * Svac(jpp, jp) - Rdag(j, jpp) * Qvac(jpp, jp);

	PSmat (j, jp) = sum;
      }
  
  // ............................
  // Calculate PS-matrix residual
  // ............................
  double Imax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ival;
	if (j == jp)
	  ival = abs (PSmat(j, jp) - complex<double> (1., 0.));
	else
	  ival = abs (PSmat(j, jp));

	if (ival > Imax)
	  Imax = ival;
      }
  
  printf ("PR, QS, and PS matrix Hermitian test residuals: %10.4e %10.4e %10.4e\n", Aamax/Ahmax, Bamax/Bhmax, Imax);

  // ...................
  // Calculate QP-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Qvac(j, jpp) * Pdag(jpp, jp);

	QPmat(j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	QPher(j, jp) = 0.5 * (QPmat(j, jp) + conj (QPmat(jp, j)));
	QPant(j, jp) = 0.5 * (QPmat(j, jp) - conj (QPmat(jp, j)));
      }
  
  // ............................
  // Calculate QP-matrix residual
  // ............................
  Ahmax = 0.; Aamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ahval = abs (QPher(j, jp));
	double aaval = abs (QPant(j, jp));

	if (ahval > Ahmax)
	  Ahmax = ahval;
	if (aaval > Aamax)
	  Aamax = aaval;	
      }

  // ...................
  // Calculate RS-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Rvac(j, jpp) * Sdag(jpp, jp);

	RSmat(j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	RSher(j, jp) = 0.5 * (RSmat(j, jp) + conj (RSmat(jp, j)));
	RSant(j, jp) = 0.5 * (RSmat(j, jp) - conj (RSmat(jp, j)));
      }
  
  // ............................
  // Calculate RS-matrix residual
  // ............................
  Bhmax = 0.; Bamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double bhval = abs (RSher(j, jp));
	double baval = abs (RSant(j, jp));

	if (bhval > Bhmax)
	  Bhmax = bhval;
	if (baval > Bamax)
	  Bamax = baval;	
      }

  // ...................
  // Calculate SP-matrix
  // ...................
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += Pvac(j, jpp) * Sdag(jpp, jp) - Qvac(j, jpp) * Rdag(jpp, jp);

	SPmat(j, jp) = sum;
      }
  
  // ............................
  // Calculate SP-matrix residual
  // ............................
  Imax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ival;
	if (j == jp)
	  ival = abs (SPmat(j, jp) - complex<double> (1., 0.));
	else
	  ival = abs (SPmat(j, jp));

	if (ival > Imax)
	  Imax = ival;
      }
  
  printf ("QP, RS, and SP matrix Hermitian test residuals: %10.4e %10.4e %10.4e\n", Aamax/Ahmax, Bamax/Bhmax, Imax);

  // ..................
  // Calculate H-matrix
  // ..................
  SolveLinearSystemTranspose (Rvac, RPmat, Pvac);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      RPdag(j, jp) = conj (RPmat(jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Hmat(j, jp) = 0.5 * RPmat(j, jp) + 0.5 * RPdag(j, jp);

  // ...................
  // Calculate iH-matrix
  // ...................
  SolveLinearSystemTranspose (Pvac, iRPmat, Rvac);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      iRPdag(j, jp) = conj (iRPmat(jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      iHmat(j, jp) = 0.5 * iRPmat(j, jp) + 0.5 * iRPdag(j, jp);
  
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

// ############################################
// Function to calculate wall boundary matrices
// ############################################
void TJ::GetVacuumWall ()
{
  printf ("Calculating wall boundary data:\n");
  
  // ...............
  // Allocate memory
  // ...............
  Pwal.resize(J, J);
  Qwal.resize(J, J);
  Rwal.resize(J, J);
  Swal.resize(J, J);

  iImat.resize(J, J);
  iIher.resize(J, J);
  iIant.resize(J, J);

  PImat.resize(J, J);
  RImat.resize(J, J);
  IRmat.resize(J, J);

  RPImat.resize(J, J);
  RPIdag.resize(J, J);
  Gmat  .resize(J, J);

  iRPImat.resize(J, J);
  iRPIdag.resize(J, J);
  iGmat.resize  (J, J);

  Rbamat.resize (J, J);

  // .........................
  // Calculate vacuum matrices
  // .........................
  int              neqns = 4*J*J;
  double           h, t_err, t;
  int              rept; rhs_chooser = 2;
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
	Pwal(j, jp) = Y[index]; index++;
	Rwal(j, jp) = Y[index]; index++;
	Qwal(j, jp) = Y[index]; index++;
	Swal(j, jp) = Y[index]; index++;
      }

  // ..........................
  // Calculate inverse I-matrix
  // ..........................
  SolveLinearSystem (Rwal, iImat, Swal);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	iIher(j, jp) = 0.5 * (iImat(j, jp) + conj (iImat(jp, j)));
	iIant(j, jp) = 0.5 * (iImat(j, jp) - conj (iImat(jp, j)));
      }

  // ....................................
  // Calculate inverse I-matrix residuals
  // ....................................
  double Ahmax = 0., Aamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double ahval = abs (iIher(j, jp));
	double aaval = abs (iIant(j, jp));

	if (ahval > Ahmax)
	  Ahmax = ahval;
	if (aaval > Aamax)
	  Aamax = aaval;	
      }
  
  printf ("I matrix Hermitian test residuals: %10.4e\n", Aamax/Ahmax);

  // .................................
  // Calculate PI, RI, and IR-matrices
  // ................................,
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum1 = - Qvac (j, jp);
	complex<double> sum2 = - Svac (j, jp);
	complex<double> sum3 = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  {
	    sum1 += Pvac(j, jpp) * iImat(jpp, jp);
	    sum2 += Rvac(j, jpp) * iImat(jpp, jp);
	    sum3 += Rvac(j, jpp) * iImat(jpp, jp);
	  }

	PImat(j, jp) = sum1;
	RImat(j, jp) = sum2;
	IRmat(j, jp) = sum3;
      }
  
  // ..................
  // Calculate G-matrix
  // ..................
  SolveLinearSystemTranspose (RImat, RPImat, PImat);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      RPIdag(j, jp) = conj (RPImat(jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Gmat(j, jp) = 0.5 * RPImat(j, jp) + 0.5 * RPIdag(j, jp);

  // ...................
  // Calculate iG-matrix
  // ...................
  SolveLinearSystemTranspose (PImat, iRPImat, RImat);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      iRPIdag(j, jp) = conj (iRPImat(jp, j));

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      iGmat(j, jp) = 0.5 * iRPImat(j, jp) + 0.5 * iRPIdag(j, jp);

  // ................
  // Calculate Rbamat
  // ................
  SolveLinearSystemTranspose (Rvac, Rbamat, Rwal);
}
 
// ####################################################
// Function to evaluate right-hand sides of vacuum odes
// ####################################################
void TJ::CashKarp45Rhs1 (double t, complex<double>* Y, complex<double>* dYdt)
{
  if (rhs_chooser == 0)
    {
      // ...........................................................
      // Right-hand sides for calculation of no-wall vacuum matrices
      // ...........................................................
      
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
  else if (rhs_chooser == 2)
    {
      // ................................................................
      // Right-hand sides for calculation of perfect-wall vacuum matrices
      // ................................................................
      
      int index = 0;
      for (int j = 0; j < J; j++)
	for (int jp = 0; jp < J; jp++)
	  {
	    int    MP  = MPOL[jp];
	    int    MMP = abs (MP);
	    
	    double m   = mpol[j];
	    double mp  = mpol[jp];
	    
	    double R    = gsl_spline_eval (Rwspline,   t, Rwacc);
	    double Z    = gsl_spline_eval (Zwspline,   t, Zwacc);
	    double R2rz = gsl_spline_eval (Rrzwspline, t, Rrzwacc);
	    double R2re = gsl_spline_eval (Rrewspline, t, Rrewacc);
	    
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


