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
  Psii.resize(J, J, NDIAG);
  Zi  .resize(J, J, NDIAG);
  Psie.resize(J, K, NDIAG);
  Ze  .resize(J, K, NDIAG);

  Umat.resize(K, K);
  Vmat.resize(K, K);
  Wmat.resize(K, K);
  Wher.resize(K, K);
  Want.resize(K, K);
  Vher.resize(K, K);
  Vant.resize(K, K);
  Uher.resize(K, K);
  Uant.resize(K, K);
  Uvec.resize(K, K);
  Ures.resize(K, K);

  Wval    = new double[K];
  Vval    = new double[K];
  Uval    = new double[K];
  deltaW  = new double[K];
  deltaWp = new double[K];
  deltaWv = new double[K];

  Psiy.resize(K, Nw+1);
  Xiy .resize(K, Nw+1);

  ya = new complex<double>[J];
  
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

  // -------------------------------------------------------------------
  // Ensure that m != 0 solutions have Z_0 = 0 at plasma-vacuum boundary
  // -------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  if (jp != jm0)
	    {
	      Psii(j, jp, i) = Psii(j, jp, i) - YYY(j,   jm0, i) * YYY(J+jm0, jp, NDIAG-1) /YYY(J+jm0, jm0, NDIAG-1);
	      Zi  (j, jp, i) = Zi  (j, jp, i) - YYY(J+j, jm0, i) * YYY(J+jm0, jp, NDIAG-1) /YYY(J+jm0, jm0, NDIAG-1);
	    }
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
  
  // ------------------------------------
  // Calculate plasma ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	Umat(j, jp) = Psii(jrnd[j], jrnd[jp], NDIAG-1);
	Vmat(j, jp) = Zi  (jrnd[j], jrnd[jp], NDIAG-1);
      }
  SolveLinearSystemTranspose (Umat, Wmat, Vmat);

  // ------------------------------------
  // Calculate vacuum ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      Vmat(j, jp) = - mpox[j] * Hmat(jrnd[j], jrnd[jp]) * mpox[jp];
  
  // -----------------------------------
  // Calculate total ideal energy matrix
  // -----------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      Umat(j, jp) = Wmat(j, jp) + Vmat(j, jp);

  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
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
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
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
  for (int j = 0; j < K; j++)
    {
      int    jmax = 0;
      double norm = -1.;

      for (int jp = 0; jp < K; jp++)
	{
	  if (abs (Uvec(jp, j)) > norm)
	    {
	      norm = abs (Uvec(jp, j));
	      jmax = jp;
	    }
	}

      complex<double> Ufac = conj (Uvec(jmax, j)) /abs (Uvec(jmax, j));

      for (int jp = 0; jp < K; jp++)
	{
	  Uvec(jp, j) *= Ufac;
	}
    }

  // ------------------------------------
  // Check orthonormality of eigenvectors
  // ------------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < K; jpp++)
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
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	double uval = abs (Ures(j, jp));

	if (uval > Umax)
	  Umax = uval;
      }

  printf ("W matrix eigenvector orthonormality test residual: %10.4e\n", Umax);

  // -----------------------------
  // Calculate ideal eigenfuctions
  // -----------------------------
  Array<complex<double>,2> Amat(K, K), Bmat(K, K);
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	Amat(j, jp) = Psii(jrnd[j], jrnd[jp], NDIAG-1);
      }

  SolveLinearSystem (Amat, Bmat, Uvec);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < K; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < K; jpp++)
	    {
	      sump += Psii(j, jrnd[jpp], i) * Bmat(jpp, jp);
	      sumz += Zi  (j, jrnd[jpp], i) * Bmat(jpp, jp);
	    }
	  
	  Psie(j, jp, i) = sump;
	  Ze  (j, jp, i) = sumz;
	}

  // ------------------------
  // Calculate delta-W values
  // ------------------------
  for (int j = 0; j < K; j++)
    {
      double sum = 0., sum1 = 0.;

      for (int jp = 0; jp < K; jp++)
	for (int jpp = 0; jpp < K; jpp++)
	  {
	    sum  += real (conj (Uvec(jp, j)) * Wmat(jp, jpp) * Uvec(jpp, j));
	    sum1 += real (conj (Uvec(jp, j)) * Vmat(jp, jpp) * Uvec(jpp, j));
	}

      deltaW [j] = M_PI*M_PI * Uval[j];
      deltaWp[j] = M_PI*M_PI * sum;
      deltaWv[j] = M_PI*M_PI * sum1;
    }

  printf ("No-wall ideal eigenvalues:\n");
  for (int j = 0; j < K; j++)
    printf ("j = %3d  deltaW = %10.3e  deltaW_p = %10.3e  deltaW_v = %10.3e\n",
    	    j, deltaW[j], deltaWp[j], deltaWv[j]);

  // --------------------------------------------------------------------------------
  // Calculate y and Z values at plasma boundary associated with ideal eigenfunctions
  // --------------------------------------------------------------------------------
  for (int j = 0; j < K; j++)
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

  // ---------------
  // Determine jzero
  // ---------------
  for (int jp = 0; jp < K; jp++)
    {
      jzero = jp;

      if (deltaWv[jp] > 0.)
	break;
    }
  
  // -------------------
  // Calculate ya values
  // -------------------
  for (int j = 0; j < J; j++)
    {
      ya[j] = Psie(j, jzero, NDIAG-1);
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
  pPsii.resize(J, J, NDIAG);
  pZi  .resize(J, J, NDIAG);
  pPsie.resize(J, K, NDIAG);
  pZe  .resize(J, K, NDIAG);

  pWmat.resize(K, K);
  pWher.resize(K, K);
  pWant.resize(K, K);
  pVmat.resize(K, K);
  pVher.resize(K, K);
  pVant.resize(K, K);
  pUmat.resize(K, K);
  pUher.resize(K, K);
  pUant.resize(K, K);
  pUvec.resize(K, K);
  pUres.resize(K, K);

  pWval    = new double[K];
  pVval    = new double[K];
  pUval    = new double[K];
  pdeltaW  = new double[K];
  pdeltaWp = new double[K];
  pdeltaWv = new double[K];

  pPsiy.resize(K, Nw+1);
  pXiy .resize(K, Nw+1);

  yb = new complex<double>[J];
  
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

  // -------------------------------------------------------------------
  // Ensure that m != 0 solutions have Z_0 = 0 at plasma-vacuum boundary
  // -------------------------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  if (jp != jm0)
	    {
	      pPsii(j, jp, i) = pPsii(j, jp, i) - YYY(j,   jm0, i) * YYY(J+jm0, jp, NDIAG-1) /YYY(J+jm0, jm0, NDIAG-1);
	      pZi  (j, jp, i) = pZi  (j, jp, i) - YYY(J+j, jm0, i) * YYY(J+jm0, jp, NDIAG-1) /YYY(J+jm0, jm0, NDIAG-1);
	    }
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
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	pUmat(j, jp) = pPsii(jrnd[j], jrnd[jp], NDIAG-1);
	pVmat(j, jp) = pZi  (jrnd[j], jrnd[jp], NDIAG-1);
      }
  SolveLinearSystemTranspose (pUmat, pWmat, pVmat);

  // ------------------------------------
  // Calculate vacuum ideal energy matrix
  // ------------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      pVmat(j, jp) = - mpox[j] * Gmat(jrnd[j], jrnd[jp]) * mpox[jp];
  
  // -----------------------------------
  // Calculate total ideal energy matrix
  // -----------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      pUmat(j, jp) = pWmat(j, jp) + pVmat(j, jp);

  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
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
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
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
  for (int j = 0; j < K; j++)
    {
      int    jmax = 0;
      double norm = -1.;

      for (int jp = 0; jp < K; jp++)
	{
	  if (abs (pUvec(jp, j)) > norm)
	    {
	      norm = abs (pUvec(jp, j));
	      jmax = jp;
	    }
	}

      complex<double> Ufac = conj (pUvec(jmax, j)) /abs (pUvec(jmax, j));

      for (int jp = 0; jp < K; jp++)
	{
	  pUvec(jp, j) *= Ufac;
	}
    }
 
  // ------------------------------------
  // Check orthonormality of eigenvectors
  // ------------------------------------
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < K; jpp++)
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
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	double uval = abs (pUres(j, jp));

	if (uval > Umax)
	  Umax = uval;
      }

  printf ("W matrix eigenvector orthonormality test residual: %10.4e\n", Umax);

  // -----------------------------
  // Calculate ideal eigenfuctions
  // -----------------------------
  Array<complex<double>,2> Amat(K, K), Bmat(K, K);
  for (int j = 0; j < K; j++)
    for (int jp = 0; jp < K; jp++)
      {
	Amat(j, jp) = pPsii(jrnd[j], jrnd[jp], NDIAG-1);
      }

  SolveLinearSystem (Amat, Bmat, pUvec);
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < K; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < K; jpp++)
	    {
	      sump += pPsii(j, jrnd[jpp], i) * Bmat(jpp, jp);
	      sumz += pZi  (j, jrnd[jpp], i) * Bmat(jpp, jp);
	    }
	  
	  pPsie(j, jp, i) = sump;
	  pZe  (j, jp, i) = sumz;
	}
 
  // ------------------------
  // Calculate delta-W values
  // ------------------------
  for (int j = 0; j < K; j++)
    {
      double sum = 0., sum1 = 0.;

      for (int jp = 0; jp < K; jp++)
	for (int jpp = 0; jpp < K; jpp++)
	  {
	    sum  += real (conj (pUvec(jp, j)) * pWmat(jp, jpp) * pUvec(jpp, j));
	    sum1 += real (conj (pUvec(jp, j)) * pVmat(jp, jpp) * pUvec(jpp, j));
	}

      pdeltaW [j] = M_PI*M_PI * pUval[j];
      pdeltaWp[j] = M_PI*M_PI * sum;
      pdeltaWv[j] = M_PI*M_PI * sum1;
    }

  printf ("Perfect-wall ideal eigenvalues:\n");
  for (int j = 0; j < K; j++)
    printf ("j = %3d  deltaW = %10.3e  deltaW_p = %10.3e  deltaW_v = %10.3e\n",
    	    j, pdeltaW[j], pdeltaWp[j], pdeltaWv[j]);

  // --------------------------------------------------------------------------------
  // Calculate y and Z values at plasma boundary associated with ideal eigenfunctions
  // --------------------------------------------------------------------------------
  for (int j = 0; j < K; j++)
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
  
  // -------------------
  // Calculate yb values
  // -------------------
  for (int j = 0; j < J; j++)
    {
      complex<double> sum = complex<double> (0., 0.);
      
      for (int jp = 0; jp < J; jp++)
	{
	  sum += Rbamat(j, jp) * mpol[jp] * ya[jp];
	}
      
      if (MPOL[j] == 0)
	yb[j] = complex<double> (0., 0.);
      else
	yb[j] = sum /mpol[j];
    }

  // ----------------
  // Determine pjzero
  // ----------------
  for (int jp = 0; jp < K; jp++)
    {
      pjzero = jp;

      if (pdeltaWv[jp] > 0.)
	break;
    }
  
  // ---------------------------
  // Calculate alphaw and gammaw
  // ---------------------------
  double sumw = 0.;
  for (int j = 0; j < J; j++)
    sumw += real (conj (yb[j]) * yb[j]);

  alphaw = M_PI*M_PI * sumw /(pdeltaWv[pjzero] - deltaWv[jzero]);
  gammaw = - deltaW[jzero] /alphaw /pdeltaW[pjzero];
  
  printf ("solution number (%1d, %1d): alphaw = %10.3e gammaw = %10.3e\n", jzero, pjzero, alphaw, gammaw);
}
