// Ideal.cpp

#include "TJ.h"

// #############################################################
// Function to calculate resistive wall mode stability of plasma
// #############################################################
void TJ::CalculateRWMStability ()
{
  // ---------------
  // Allocate memory
  // ---------------
  Wnw.resize (J, J);
  Wpw.resize (J, J);
  Wpwh.resize(J, J);
  Wpw2.resize(J, J);

  Dmat.resize(J, J);
  Dher.resize(J, J);
  Dsqt.resize(J, J);
  Dinv.resize(J, J);
  
  Ehmt.resize(J, J);

  FFmt.resize(J, J);
  FFhr.resize(J, J);
  FFan.resize(J, J);

  FFvc.resize(J, J);
  FFrs.resize(J, J);

  FFvl = new double[J];
		   
  // --------------------------------
  // Calculate W_nw and W_pw matrices
  // --------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Wnw(j, jp) = Wmat(j, jp) - Hmat(j, jp);
	Wpw(j, jp) = Wmat(j, jp) - Gmat(j, jp);
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Wpwh(j, jp) = 0.5 * (Wpw(j, jp) + conj (Wpw(jp, j)));
      }

  // ------------------
  // Calculate D-matrix
  // ------------------
  InvSquareRootMatrix (Wpwh, Wpw2);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp  = 0; jpp < J; jpp++)
	  for (int jppp = 0; jppp < J; jppp++)
	    sum += Wpw2(j, jpp) * Cher(jpp, jppp) * Wpw2(jppp, jp);

	Dmat(j, jp) = sum;
      }
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Dher(j, jp) = 0.5 * (Dmat (j, jp) + conj (Dmat (jp, j)));
      }

  // ------------------
  // Calculate E-matrix
  // ------------------
  SquareRootMatrix    (Dher, Dsqt);
  InvSquareRootMatrix (Dher, Dinv);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);
	
	for (int jpp = 0; jpp < J; jpp++)
	  for (int jppp = 0; jppp < J; jppp++)
	    for (int jpppp = 0; jpppp < J; jpppp++)
	      sum += Wpw2(j, jpp) * Wnw(jpp, jppp) * Bmat(jppp, jpppp) * Wpw2(jpppp, jp);
	
	Ehmt(j, jp) = sum;
      }

  // ------------------
  // Calculate F-matrix
  // ------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  for (int jppp = 0; jppp < J; jppp++)
	    sum += Dsqt(j, jpp) * Ehmt(jpp, jppp) * Dinv(jppp, jp);

	FFmt(j, jp) = sum;
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	FFhr(j, jp) = 0.5 * (FFmt (j, jp) + conj (FFmt (jp, j)));
	FFan(j, jp) = 0.5 * (FFmt (j, jp) - conj (FFmt (jp, j)));
      }

  // ...........................
  // Calculate F-matrix residual
  // ...........................
  double Fhmax = 0., Famax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double fhval = abs (FFhr (j, jp));
	double faval = abs (FFan (j, jp));

	if (fhval > Fhmax)
	  Fhmax = fhval;
	if (faval > Famax)
	  Famax = faval;	
      }

  printf ("F matrix Hermitian test residual: %10.4e\n", Famax/Fhmax);

  // --------------------------------------------------
  // Calculate eigenvalues and eigenvectors of F-matrix
  // --------------------------------------------------
  GetEigenvalues (FFhr, FFvl, FFvc);
  
  // --------------------------------------------------------------------
  // Adjust eigenvectors such that element with largest magnitude is real
  // --------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    {
      int    jmax = 0;
      double norm = -1.;

      for (int jp = 0; jp < J; jp++)
	{
	  if (abs (FFvc(jp, j)) > norm)
	    {
	      norm = abs (FFvc(jp, j));
	      jmax = jp;
	    }
	}

      complex<double> Ffac = conj (FFvc(jmax, j)) /abs (FFvc(jmax, j));

      for (int jp = 0; jp < J; jp++)
	{
	  FFvc(jp, j) *= Ffac;
	}
    }

  // ------------------------------------
  // Check orthonormality of eigenvectors
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += conj (FFvc(jpp, j)) * FFvc(jpp, jp);

	if (j == jp)
	  FFrs(j, jp) = sum - complex<double> (1., 0.);
	else
	  FFrs(j, jp) = sum;
      }

  // ----------------------------------
  // Calculate orthonormality residuals
  // ----------------------------------
  double Fmax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double fval = abs (FFrs(j, jp));

	if (fval > Fmax)
	  Fmax = fval;
      }

  printf ("F matrix eigenvector orthonormality test residual: %10.4e\n", Fmax);

  // ---------------------------------------
  // Display resistive wall mode eigenvalues
  // ---------------------------------------
  printf ("Resistive wall mode eigenvalues:\n");
  for (int j = 0; j < J; j++)
    printf ("j = %3d  fw = %10.3e\n", j, -FFvl[j]);
}

