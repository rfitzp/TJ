// Ideal.cpp

#include "TJ.h"

// ###############################################
// Function to calculate ideal stabiltiy of plasma
// ###############################################
void TJ::CalculateIdealStability ()
{
  // ---------------
  // Allocate memory
  // ---------------
  Psii .resize(J, J, NDIAG);
  Zi   .resize(J, J, NDIAG);
  Xii  .resize(J, J, NDIAG);
  Chii .resize(J, J, NDIAG);
  Ji   .resize(J, J);
  Wmat .resize(J, J);
  Wher .resize(J, J);
  Want .resize(J, J);
  Vmat .resize(J, J);
  Vher .resize(J, J);
  Vant .resize(J, J);
  Umat .resize(J, J);
  Uher .resize(J, J);
  Uant .resize(J, J);
  Uvec .resize(J, J);
  Ures .resize(J, J);
  Psie .resize(J, J, NDIAG);
  Ze   .resize(J, J, NDIAG);
  Xie  .resize(J, J, NDIAG);
  Je   .resize(J, J);
  lvals.resize(J, NDIAG);

  Wval    = new double[J];
  Vval    = new double[J];
  Uval    = new double[J];
  deltaW  = new double[J];
  deltaWp = new double[J];
  deltaWv = new double[J];

  Psiy.resize(J, Nw+1);
  Jy  .resize(J, Nw+1);
  Xiy .resize(J, Nw+1);

  gamma  = new complex<double>[J];
  gammax = new complex<double>[J];
  
  // -----------------------------------------------------
  // Calculate ideal solutions launched from magnetic axis
  // -----------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = YYY(j,   jp, i);
	  complex<double> sumz = YYY(J+j, jp, i);
	  double          km   = Getkm (Rgrid[i], MPOL[j]);
	  
	  for (int k = 0; k < nres; k++)
	    {
	      sump -= Psiu(j, k, i) * Pia(k, jp);
	      sumz -= Zu  (j, k, i) * Pia(k, jp);
	    }
	  
	  Psii(j, jp, i) = sump;
	  Zi  (j, jp, i) = sumz;
	  Xii (j, jp, i) = sump /(mpol[j] - ntor * Getq (Rgrid[i]));
	  Chii(j, jp, i) = km * sump + sumz;
	}
  
  // -----------------------------------------------------
  // Normalize ideal solutions launched from magnetic axis
  // -----------------------------------------------------
  for (int jp = 0; jp < J; jp++)
    {
      double norm = 0.;

      for (int j = 0; j < J; j++)
	norm += real (conj (Psii(j, jp, NDIAG-1)) * Psii(j, jp, NDIAG-1));

      for (int i = 0; i < NDIAG; i++)
	for (int j = 0; j < J; j++)
	  {
	    Psii(j, jp, i) /= sqrt(norm);
	    Zi  (j, jp, i) /= sqrt(norm);
	    Xii (j, jp, i) /= sqrt(norm);
	    Chii(j, jp, i) /= sqrt(norm);
	  }
    }
  
  // ---------------------------------------------------------------------------------------
  // Calculate boundary currents associated with ideal solutions launched from magnetic axis
  // ---------------------------------------------------------------------------------------
  double qa = Getq (1.);
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = Zi(j, jp, NDIAG-1) /(mpol[j] - ntor*qa);

	for (int jpp = 0; jpp < J; jpp++)
	  if (FREE > 0 || FREE < 0)
	    sum -= Hmat(j, jpp) * Psii(jpp, jp, NDIAG-1);
	  else
	    sum -= Gmat(j, jpp) * Psii(jpp, jp, NDIAG-1);

	Ji(j, jp) = sum;
      }

  // ------------------------------------
  // Calculate plasma ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Umat(j, jp) = Psii(j, jp, NDIAG-1);
	Vmat(j, jp) = Zi(j, jp, NDIAG-1) /(mpol[j] - ntor*qa);
      }
  SolveLinearSystemTranspose (Umat, Wmat, Vmat);

  // ------------------------------------
  // Calculate vacuum ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      if (FREE > 0 || FREE < 0)
	Vmat(j, jp) = - Hmat(j, jp);
      else
	Vmat(j, jp) = - Gmat(j, jp);

  // -----------------------------------
  // Calculate total ideal energy matrix
  // -----------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Umat(j, jp) = Wmat(j, jp) + Vmat(j, jp);

  // -------------------------------
  // Transform matrices if XI is set
  // -------------------------------
  if (XI)
    for (int j = 0; j < J; j++)
      for (int jp = 0; jp < J; jp++)
	{
	  Wmat(j, jp) = (mpol[j] - ntor*qa) * Wmat(j, jp) * (mpol[jp] - ntor*qa);
	  Vmat(j, jp) = (mpol[j] - ntor*qa) * Vmat(j, jp) * (mpol[jp] - ntor*qa);
	  Umat(j, jp) = (mpol[j] - ntor*qa) * Umat(j, jp) * (mpol[jp] - ntor*qa);
	}
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
 	Wher(j, jp) = 0.5 * (Wmat(j, jp) + conj (Wmat(jp, j)));
	Want(j, jp) = 0.5 * (Wmat(j, jp) - conj (Wmat(jp, j)));
	Vher(j, jp) = 0.5 * (Vmat(j, jp) + conj (Vmat(jp, j)));
	Vant(j, jp) = 0.5 * (Vmat(j, jp) - conj (Vmat(jp, j)));
	Uher(j, jp) = 0.5 * (Umat(j, jp) + conj (Umat(jp, j)));
	Uant(j, jp) = 0.5 * (Umat(j, jp) - conj (Umat(jp, j)));
      }

  // --------------------------
  // Calculate matrix residuals
  // --------------------------
  double Whmax = 0., Wamax = 0.;
  double Vhmax = 0., Vamax = 0.;
  double Uhmax = 0., Uamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double whval = abs (Wher(j, jp));
	double waval = abs (Want(j, jp));
	double vhval = abs (Vher(j, jp));
	double vaval = abs (Vant(j, jp));
	double uhval = abs (Uher(j, jp));
	double uaval = abs (Uant(j, jp));

	if (whval > Whmax)
	  Whmax = whval;
	if (waval > Wamax)
	  Wamax = waval;
	if (vhval > Vhmax)
	  Vhmax = vhval;
	if (vaval > Vamax)
	  Vamax = vaval;
	if (uhval > Uhmax)
	  Uhmax = uhval;
	if (uaval > Uamax)
	  Uamax = uaval;
      }

  printf ("Wp, Wv, and W matrix Hermitian test residuals: %10.4e %10.4e %10.4e\n", Wamax /Whmax, Vamax /Vhmax, Uamax /Uhmax);
 
  // ------------------------------------------
  // Calculate eigenvalues of W- and V-matrices
  // ------------------------------------------
  GetEigenvalues (Wher, Wval);
  GetEigenvalues (Vher, Vval);
  
  // --------------------------------------------------
  // Calculate eigenvalues and eigenvectors of U-matrix
  // --------------------------------------------------
  GetEigenvalues (Uher, Uval, Uvec);
  
  // --------------------------------------------------------------------
  // Adjust eigenvectors such that element with largest magnitude is real
  // --------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    {
      int    jmax = 0;
      double norm = -1.;

      for (int jp = 0; jp < J; jp++)
	{
	  if (abs (Uvec(jp, j)) > norm)
	    {
	      norm = abs (Uvec(jp, j));
	      jmax = jp;
	    }
	}

      complex<double> Ufac = conj (Uvec(jmax, j)) /abs (Uvec(jmax, j));

      for (int jp = 0; jp < J; jp++)
	{
	  Uvec(jp, j) *= Ufac;
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
	  sum += conj (Uvec(jpp, j)) * Uvec(jpp, jp);

	if (j == jp)
	  Ures(j, jp) = sum - complex<double> (1., 0.);
	else
	  Ures(j, jp) = sum;
      }

  // ----------------------------------
  // Calculate orthonormality residuals
  // ----------------------------------
  double Umax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double uval = abs (Ures(j, jp));

	if (uval > Umax)
	  Umax = uval;
      }

  printf ("W matrix eigenvector orthonormality test residual: %10.4e\n", Umax);

  // -----------------------------
  // Calculate ideal eigenfuctions
  // -----------------------------
  Array<complex<double>,2> Amat(J, J), Bmat(J, J);
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	if (XI)
	  Amat(j, jp) = Xii(j, jp, NDIAG-1);
	else
	  Amat(j, jp) = Psii(j, jp, NDIAG-1);
      }

  SolveLinearSystem (Amat, Bmat, Uvec);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  complex<double> sumx = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < J; jpp++)
	    {
	      sump += Psii(j, jpp, i) * Bmat(jpp, jp);
	      sumz += Zi  (j, jpp, i) * Bmat(jpp, jp);
	      sumx += Xii (j, jpp, i) * Bmat(jpp, jp);
	    }
	  
	  Psie(j, jp, i) = sump;
	  Ze  (j, jp, i) = sumz;
	  Xie (j, jp, i) = sumx;
	}

  // ----------------------------------------------------------------
  // Calculate boundary currents associated with ideal eigenfunctions
  // ----------------------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	if (XI)
	  Je(j, jp) = Uvec(j, jp) * Uval[jp] /(mpol[j] - ntor*qa);
	else
	  Je(j, jp) = Uvec(j, jp) * Uval[jp];
      }

  // ------------------------
  // Calculate delta-W values
  // ------------------------
  for (int j = 0; j < J; j++)
    {
      double sum = 0., sum1 = 0.;;

      for (int jp = 0; jp < J; jp++)
	for (int jpp = 0; jpp < J; jpp++)
	  {
	    sum  += real (conj (Uvec(jp, j)) * Wmat(jp, jpp) * Uvec(jpp, j));
	    sum1 += real (conj (Uvec(jp, j)) * Vmat(jp, jpp) * Uvec(jpp, j));
	}

      deltaW [j] = M_PI*M_PI * Uval[j];
      deltaWp[j] = M_PI*M_PI * sum;
      deltaWv[j] = M_PI*M_PI * sum1;
    }

  printf ("Ideal eigenvalues:\n");
  for (int j = 0; j < J; j++)
    printf ("j = %3d  deltaW = %10.3e  deltaW_p = %10.3e  deltaW_v = %10.3e\n",
    	    j, deltaW[j], deltaWp[j], deltaWv[j]);

  // ----------------------------------------------------------------------------------
  // Calculate Psi and J values at plasma boundary associated with ideal eigenfunctions
  // ----------------------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    {
      for (int i = 0; i <= Nw; i++)
	{
	  double theta = tbound[i];

	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumj = complex<double> (0., 0.);
	  complex<double> sumx = complex<double> (0., 0.);

	  for (int jp = 0; jp < J; jp++)
	    {
	      sump += Psie(jp, j, NDIAG-1)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	      sumj += Je(jp, j)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	      sumx += Xie(jp, j, NDIAG-1)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	    }

	  Psiy(j, i) = sump;
	  Jy  (j, i) = sumj;
	  Xiy (j, i) = sumx;
	}
    }

  // -------------------------------------------------------------------------------------
  // Calculate expansion of Psi_x and Psi_rmp at boundary in terms of ideal eigenfunctions
  // -------------------------------------------------------------------------------------
  if (RMP)
    {
      for (int j = 0; j < J; j++)
	{
	  complex<double> sumx = complex<double> (0., 0.);
	  complex<double> sum  = complex<double> (0., 0.);
	  
	  for (int jp = 0; jp < J; jp++)
	    {
	      if (XI)
		{
		  sumx += conj (Uvec(jp, j)) * Psix[jp]            /(mpol[jp] - ntor*qa);
		  sum  += conj (Uvec(jp, j)) * Psirmp(jp, NDIAG-1) /(mpol[jp] - ntor*qa);
		}
	      else
		{
		  sumx += conj (Uvec(jp, j)) * Psix[jp];
		  sum  += conj (Uvec(jp, j)) * Psirmp(jp, NDIAG-1);
		}
	    }
	  
	  gammax[j] = sumx;
	  gamma [j] = sum;
	}
    }
  
  // ------------------------------------------------------
  // Calculate eigenvalues of plasma energy matrix versus r
  // ------------------------------------------------------
  if (INTR)
    {
      for (int i = 1; i < NDIAG; i++)
	{
	  Array<complex<double>,2> Xiimat(J, J);
	  Array<complex<double>,2> Chimat(J, J);
	  Array<complex<double>,2> Emat  (J, J);
	  Array<complex<double>,2> hEmat (J, J);
	  
	  double* evals = new double[J];
	  
	  for (int j = 0; j < J; j++)
	    for (int jp = 0; jp < J; jp++)
	      {
		Xiimat(j, jp) = Xii (j, jp, i);
		Chimat(j, jp) = Chii(j, jp, i);
	      }
	  
	  SolveLinearSystemTranspose (Xiimat, Emat, Chimat);
	  
	  for (int j = 0; j < J; j++)
	    for (int jp = 0; jp < J; jp++)
	      {
		hEmat(j, jp) = (Emat(j, jp) + conj (Emat(jp, j))) /2.;
	      }
	  
	  GetEigenvalues (hEmat, evals);
	  
	  for (int j = 0; j < J; j++)
	    lvals(j, i) = M_PI*M_PI * evals[j];
	  
	  delete[] evals;
	}
      
      for (int j = 0; j < J; j++)
	lvals(j, 0) = lvals(j, 1);
    }
}
