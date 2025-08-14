// Ideal.cpp

#include "Vertical.h"

// #######################################################
// Function to calculate no-wall ideal stability of plasma
// #######################################################
void Vertical::CalculateNoWallIdealStability ()
{
  // ---------------
  // Allocate memory
  // ---------------
  Psii .resize(J, J, NDIAG);
  Zi   .resize(J, J, NDIAG);
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

  Wval    = new double[J];
  Vval    = new double[J];
  Uval    = new double[J];
  deltaW  = new double[J];
  deltaWp = new double[J];
  deltaWv = new double[J];

  Psiy.resize(J, Nw+1);
  Xiy .resize(J, Nw+1);

  // -----------------------------------------
  // Get solutions launched from magnetic axis
  // -----------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  Psii(j, jp, i) = YYY(j,   jp, i);
	  Zi  (j, jp, i) = YYY(J+j, jp, i);
	}

  // -----------------------------------------------
  // Normalize solutions launched from magnetic axis
  // -----------------------------------------------
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
	  }
    }
  
  // -----------------------------------
  // Calculate plasma ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Umat(j, jp) = Psii(j, jp, NDIAG-1);
	Vmat(j, jp) = Zi  (j, jp, NDIAG-1);
      }
  SolveLinearSystemTranspose (Umat, Wmat, Vmat);

  // ------------------------------------
  // Calculate vacuum ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Vmat(j, jp) = - mpol[j] * Hmat(j, jp) * mpol[jp];
  
  // -----------------------------------
  // Calculate total ideal energy matrix
  // -----------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      Umat(j, jp) = Wmat(j, jp) + Vmat(j, jp);

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
	Amat(j, jp) = Psii(j, jp, NDIAG-1);
      }

  SolveLinearSystem (Amat, Bmat, Uvec);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < J; jpp++)
	    {
	      sump += Psii(j, jpp, i) * Bmat(jpp, jp);
	      sumz += Zi  (j, jpp, i) * Bmat(jpp, jp);
	    }
	  
	  Psie(j, jp, i) = sump;
	  Ze  (j, jp, i) = sumz;
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

  printf ("No-wall ideal eigenvalues:\n");
  for (int j = 0; j < J; j++)
    printf ("j = %3d  deltaW = %10.3e  deltaW_p = %10.3e  deltaW_v = %10.3e\n",
    	    j, deltaW[j], deltaWp[j], deltaWv[j]);

  // --------------------------------------------------------------------------------
  // Calculate y and Z values at plasma boundary associated with ideal eigenfunctions
  // --------------------------------------------------------------------------------
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
	      sumx += Ze(jp, j, NDIAG-1)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	    }

	  Psiy(j, i) = sump;
	  Xiy (j, i) = sumx;
	}
    }
}

// ############################################################
// Function to calculate perfect-wall ideal stability of plasma
// ############################################################
void Vertical::CalculatePerfectWallIdealStability ()
{
  // ---------------
  // Allocate memory
  // ---------------
  pPsii .resize(J, J, NDIAG);
  pZi   .resize(J, J, NDIAG);
  pWmat .resize(J, J);
  pWher .resize(J, J);
  pWant .resize(J, J);
  pVmat .resize(J, J);
  pVher .resize(J, J);
  pVant .resize(J, J);
  pUmat .resize(J, J);
  pUher .resize(J, J);
  pUant .resize(J, J);
  pUvec .resize(J, J);
  pUres .resize(J, J);
  pPsie .resize(J, J, NDIAG);
  pZe   .resize(J, J, NDIAG);

  pWval    = new double[J];
  pVval    = new double[J];
  pUval    = new double[J];
  pdeltaW  = new double[J];
  pdeltaWp = new double[J];
  pdeltaWv = new double[J];

  pPsiy.resize(J, Nw+1);
  pXiy .resize(J, Nw+1);

  // -----------------------------------------
  // Get solutions launched from magnetic axis
  // -----------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  pPsii(j, jp, i) = YYY(j,   jp, i);
	  pZi  (j, jp, i) = YYY(J+j, jp, i);
	}

  // -----------------------------------------------
  // Normalize solutions launched from magnetic axis
  // -----------------------------------------------
  for (int jp = 0; jp < J; jp++)
    {
      double norm = 0.;

      for (int j = 0; j < J; j++)
	norm += real (conj (pPsii(j, jp, NDIAG-1)) * pPsii(j, jp, NDIAG-1));

      for (int i = 0; i < NDIAG; i++)
	for (int j = 0; j < J; j++)
	  {
	    pPsii(j, jp, i) /= sqrt(norm);
	    pZi  (j, jp, i) /= sqrt(norm);
	  }
    }
  
  // -----------------------------------
  // Calculate plasma ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	pUmat(j, jp) = pPsii(j, jp, NDIAG-1);
	pVmat(j, jp) = pZi  (j, jp, NDIAG-1);
      }
  SolveLinearSystemTranspose (pUmat, pWmat, pVmat);

  // ------------------------------------
  // Calculate vacuum ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      pVmat(j, jp) = - mpol[j] * Gmat(j, jp) * mpol[jp];
  
  // -----------------------------------
  // Calculate total ideal energy matrix
  // -----------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      pUmat(j, jp) = pWmat(j, jp) + pVmat(j, jp);

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
 	pWher(j, jp) = 0.5 * (pWmat(j, jp) + conj (pWmat(jp, j)));
	pWant(j, jp) = 0.5 * (pWmat(j, jp) - conj (pWmat(jp, j)));
	pVher(j, jp) = 0.5 * (pVmat(j, jp) + conj (pVmat(jp, j)));
	pVant(j, jp) = 0.5 * (pVmat(j, jp) - conj (pVmat(jp, j)));
	pUher(j, jp) = 0.5 * (pUmat(j, jp) + conj (pUmat(jp, j)));
	pUant(j, jp) = 0.5 * (pUmat(j, jp) - conj (pUmat(jp, j)));
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
	double whval = abs (pWher(j, jp));
	double waval = abs (pWant(j, jp));
	double vhval = abs (pVher(j, jp));
	double vaval = abs (pVant(j, jp));
	double uhval = abs (pUher(j, jp));
	double uaval = abs (pUant(j, jp));

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
  GetEigenvalues (pWher, pWval);
  GetEigenvalues (pVher, pVval);
  
  // --------------------------------------------------
  // Calculate eigenvalues and eigenvectors of U-matrix
  // --------------------------------------------------
  GetEigenvalues (pUher, pUval, pUvec);
  
  // --------------------------------------------------------------------
  // Adjust eigenvectors such that element with largest magnitude is real
  // --------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    {
      int    jmax = 0;
      double norm = -1.;

      for (int jp = 0; jp < J; jp++)
	{
	  if (abs (pUvec(jp, j)) > norm)
	    {
	      norm = abs (pUvec(jp, j));
	      jmax = jp;
	    }
	}

      complex<double> Ufac = conj (pUvec(jmax, j)) /abs (pUvec(jmax, j));

      for (int jp = 0; jp < J; jp++)
	{
	  pUvec(jp, j) *= Ufac;
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
	  sum += conj (pUvec(jpp, j)) * pUvec(jpp, jp);

	if (j == jp)
	  pUres(j, jp) = sum - complex<double> (1., 0.);
	else
	  pUres(j, jp) = sum;
      }

  // ----------------------------------
  // Calculate orthonormality residuals
  // ----------------------------------
  double Umax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double uval = abs (pUres(j, jp));

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
	Amat(j, jp) = pPsii(j, jp, NDIAG-1);
      }

  SolveLinearSystem (Amat, Bmat, pUvec);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < J; jpp++)
	    {
	      sump += pPsii(j, jpp, i) * Bmat(jpp, jp);
	      sumz += pZi  (j, jpp, i) * Bmat(jpp, jp);
	    }
	  
	  pPsie(j, jp, i) = sump;
	  pZe  (j, jp, i) = sumz;
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
	    sum  += real (conj (pUvec(jp, j)) * pWmat(jp, jpp) * pUvec(jpp, j));
	    sum1 += real (conj (pUvec(jp, j)) * pVmat(jp, jpp) * pUvec(jpp, j));
	}

      pdeltaW [j] = M_PI*M_PI * pUval[j];
      pdeltaWp[j] = M_PI*M_PI * sum;
      pdeltaWv[j] = M_PI*M_PI * sum1;
    }

  printf ("Perfect-wall ideal eigenvalues:\n");
  for (int j = 0; j < J; j++)
    printf ("j = %3d  deltaW = %10.3e  deltaW_p = %10.3e  deltaW_v = %10.3e\n",
    	    j, pdeltaW[j], pdeltaWp[j], pdeltaWv[j]);

  // --------------------------------------------------------------------------------
  // Calculate y and Z values at plasma boundary associated with ideal eigenfunctions
  // --------------------------------------------------------------------------------
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
	      sump += pPsie(jp, j, NDIAG-1)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	      sumx += pZe(jp, j, NDIAG-1)
		* complex<double> (cos (mpol[jp] * theta), sin (mpol[jp] * theta));
	    }

	  pPsiy(j, i) = sump;
	  pXiy (j, i) = sumx;
	}
    }
}
