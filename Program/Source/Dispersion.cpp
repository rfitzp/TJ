// Dispersion.cpp

#include "TJ.h"

// ######################################################################################
// Function to find tearing mode dispersion relation and construct tearing eigenfunctions
// ######################################################################################
void TJ::FindDispersion ()
{
  printf ("Dispersion relation data:\n");

  // .............
  // Assign memory
  // .............
  Psia .resize(J,    J);
  Za   .resize(J,    J);
  Pia  .resize(nres, J);
  Psis .resize(J,    nres);
  Zs   .resize(J,    nres);
  Pis  .resize(nres, nres);
  Xmat .resize(J,    J);
  Ymat .resize(J,    nres);
  Omat .resize(J,    nres);
  Fmat .resize(nres, nres);
  Fher .resize(nres, nres);
  Fvec .resize(nres, nres);
  Emat .resize(nres, nres);
  Psif .resize(J,    nres, NDIAG);
  Zf   .resize(J,    nres, NDIAG);
  Tf   .resize(nres, NDIAG);
  Psiu .resize(J,    nres, NDIAG);
  Zu   .resize(J,    nres, NDIAG);
  Tu   .resize(nres, NDIAG);
  Tfull.resize(nres, nres, NDIAG);
  Tunrc.resize(nres, nres, NDIAG);
  Fval = new double[nres];

  // ............
  // Collate data
  // ............
  int index = NDIAG - 1;
  
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Psia(j, jp) = YYY(j,   jp, index);
	Za  (j, jp) = YYY(J+j, jp, index);
      }

  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Pia(j, jp) = Pi(j, jp);
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Psis(j, jp) = YYY(j,   J+jp, index);
	Zs  (j, jp) = YYY(J+j, J+jp, index);
      }

  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Pis(j, jp) = Pi(j, J+jp);
      }

  // ...........................
  // Construct X- and Y-matrices
  // ...........................
  double qa = Getq (1.);
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	if (FREE)
	  Xmat(j, jp) = Za(j, jp) /(mpol[j] - ntor*qa);
	else
	  Xmat(j, jp) = complex<double> (0., 0.);

	for (int k = 0; k < J; k++)
	  Xmat(j, jp) -= Hmat(j, k) * Psia(k, jp);
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	if (FREE)
	  Ymat(j, jp) = - Zs(j, jp) /(mpol[j] - ntor*qa) /dPi(jp);
	else
	  Ymat(j, jp) = complex<double> (0., 0.);

	for (int k = 0; k < J; k++)
	  Ymat(j, jp) += Hmat(j, k) * Psis(k, jp) /dPi(jp);
      }

  // ......................
  // Calculate Omega-matrix
  // ......................
  SolveLinearSystem (Xmat, Omat, Ymat);

  // ..................
  // Calculate F-matrix
  // ..................
  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	Fmat(j, jp) = Pis(j, jp) /dPi(jp);

	for (int k = 0; k < J; k++)
	  Fmat(j, jp) += Pia(j, k) * Omat(k, jp);
      }

  // ...............
  // Output F-matrix
  // ...............
  printf ("Re(F):\n");
  for (int j = 0; j < nres; j++)
    {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", real(Fmat(j, jp)));
      printf ("\n");
    }
  printf ("Im(F):\n");
  for (int j = 0; j < nres; j++)
   {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", imag(Fmat(j, jp)));
      printf ("\n");
    }

  // ..................................................
  // Calculate eigenvalues and eigenvectors of F-matrix
  // ..................................................
  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      Fher (j, jp) = 0.5 * (Fmat (j, jp) + conj (Fmat (jp, j)));

  GetEigenvalues (Fher, Fval, Fvec);

  printf ("Eigenvalues of F-matrix:");
  for (int i = 0; i < nres; i++)
    printf (" %10.3e ", Fval[i]);
  printf ("\n");

  printf ("Eigenvectors of F-matrix:\n");
  for (int i = 0; i < nres; i++)
    {
      for (int j = 0; j < nres; j++)
	printf ("(%10.3e, %10.3e) ", real (Fvec(i, j)), imag (Fvec(i, j)));
      printf ("\n");
    }
  
  // ..................................................
  // Calculate fully-reconnected tearing eigenfunctions
  // ..................................................
  for (int i = 0; i < NDIAG; i++)
    {
      for (int j = 0; j < J; j++)
	for (int k = 0; k < nres; k++)
	  {
	    Psif(j, k, i) = YYY(j, J+k, i) /dPi(k);

	    for (int jp = 0; jp < J; jp++)
	      Psif(j, k, i) += YYY(j, jp, i) * Omat(jp, k);

	    Zf(j, k, i) = YYY(J+j, J+k, i) /dPi(k);

	    for (int jp = 0; jp < J; jp++)
	      Zf(j, k, i) += YYY(J+j, jp, i) * Omat(jp, k);
	  }
    }

  // ..................................................................
  // Calculate torques associated with fully-reconnected eigenfunctions
  // ..................................................................
  for (int i = 0; i < NDIAG; i++)
    for (int k = 0; k < nres; k++)
      {
	complex<double> I   = complex<double> (0., 1.);
	complex<double> Sum = complex<double> (0., 0.);
	double          q   = Getq (Rgrid[i]);
	
	for (int j = 0; j < J; j++)
	  {
	    double mj  = mpol[j];
	    double mnq = mj - ntor * q;
	    
	    Sum += (conj(Zf(j, k, i)) * Psif(j, k, i) - conj(Psif(j, k, i)) * Zf(j, k, i)) /mnq;
	  }
	Sum *= I * M_PI*M_PI * ntor;

	Tf(k, i) = real(Sum);
      }

  // ...........................................................................
  // Calculate torques associated with pairs of fully-reconnected eigenfunctions
  // ...........................................................................
  GetTorqueFull ();

  // ..................
  // Calculate E-matrix
  // ..................
  InvertMatrix (Fmat, Emat);

  // ...............
  // Output E-matrix
  // ...............
  printf ("Re(E):\n");
  for (int j = 0; j < nres; j++)
    {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", real (Emat(j, jp)));
      printf ("\n");
    }
  printf ("Im(E):\n");
  for (int j = 0; j < nres; j++)
   {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", imag (Emat(j, jp)));
      printf ("\n");
    } 
  printf ("Re(E_res):\n");
  for (int j = 0; j < nres; j++)
    {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", real (Emat(j, jp) - conj (Emat(jp, j))));
      printf ("\n");
    }
  printf ("Im(E_res):\n");
  for (int j = 0; j < nres; j++)
   {
      for (int jp = 0; jp < nres; jp++)
	printf ("%10.3e ", imag (Emat(j, jp) - conj (Emat(jp, j))));
      printf ("\n");
    }

  double Eerr = 0.;
  for (int j = 0; j < nres; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	double eval = abs (Emat(j, jp) - conj (Emat(jp, j)));

	if (eval > Eerr)
	  Eerr = eval;	
      }
  printf ("E-matrix Hermitian test: %10.3e\n", Eerr);

  // ..............................................
  // Calculate unreconnected tearing eigenfunctions
  // ..............................................
  for (int i = 0; i < NDIAG; i++)
    {
      for (int j = 0; j < J; j++)
	for (int k = 0; k < nres; k++)
	  {
	    Psiu(j, k, i) = complex<double> (0., 0.);

	    for (int kp = 0; kp < nres; kp++)
	      Psiu(j, k, i) += Psif(j, kp, i) * Emat(kp, k);

	    Zu(j, k, i) = complex<double> (0., 0.);

	    for (int kp = 0; kp < nres; kp++)
	      Zu(j, k, i) += Zf(j, kp, i) * Emat(kp, k);
	  }
    }

  // ..............................................................
  // Calculate torques associated with unreconnected eigenfunctions
  // ..............................................................
  for (int i = 0; i < NDIAG; i++)
    for (int k = 0; k < nres; k++)
      {
	complex<double> I   = complex<double> (0., 1.);
	complex<double> Sum = complex<double> (0., 0.);
	double          q   = Getq (Rgrid[i]);
	
	for (int j = 0; j < J; j++)
	  {
	    double mj  = mpol[j];
	    double mnq = mj - ntor * q;
	    
	    Sum += (conj(Zu(j, k, i)) * Psiu(j, k, i) - conj(Psiu(j, k, i)) * Zu(j, k, i)) /mnq;
	  }
	Sum *= I * M_PI*M_PI * ntor;

	Tu(k, i) = real(Sum);
      }

  // .......................................................................
  // Calculate torques associated with pairs of unreconnected eigenfunctions
  // .......................................................................
  GetTorqueUnrc ();
}

// #########################################################
// Function to calculate resonant magnetic perturbation data
// #########################################################
void TJ::CalculateResonantMagneticPerturbation ()
{
  // ...............
  // Allocate memory
  // .............
  Upsilon = new complex<double>[J];
  Lambda  = new complex<double>[J];
  Chi     = new complex<double>[J];

  Psirmp.resize(J, NDIAG);
  Zrmp  .resize(J, NDIAG);

  // ........................
  // Calculate Upsilon-vector
  // ........................
  SolveLinearSystem (Xmat, Upsilon, Xi);

  // .......................
  // Calculate Lambda-vector
  // .......................
  for (int k = 0; k < nres; k++)
    {
      complex<double> sum = complex<double> (0., 0.);

      for (int j = 0; j < J; j++)
	sum += Pia(k, j) * Upsilon[j];

      Lambda[k] = sum;
    }

  // ....................
  // Calculate Chi-vector
  // ....................
  for (int k = 0; k < nres; k++)
    {
      complex<double> sum = complex<double> (0., 0.);

      for (int kp = 0; kp < nres; kp++)
	sum += Emat(k, kp) * Lambda[kp];
      
      Chi[k] = sum;
    }

  // .................
  // Output Chi-vector
  // .................
  printf ("Chi vector:\n");
  for (int k = 0; k < nres; k++)
    printf ("Rational surface %2d: Chi = (%10.3e, %10.3e) |Chi| = %10.3e\n",
	    k+1, real(Chi[k]), imag(Chi[k]), abs(Chi[k]));

  // .........................
  // Calculate Psirmp and Zrmp
  // .........................
  for (int i = 0; i < NDIAG; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jp = 0; jp < J; jp++)
	    {
	      complex<double> sump1 = - YYY(j,   jp, i);
	      complex<double> sumz1 = - YYY(J+j, jp, i);
	      
	      for (int k = 0; k < nres; k++)
		{
		  sump1 += Psiu(j, k, i) * Pia(k, jp);
		  sumz1 += Zu  (j, k, i) * Pia(k, jp);
		}

	      sump += sump1 * Upsilon[jp];
	      sumz += sumz1 * Upsilon[jp];
	    }

	  Psirmp(j, i) = sump;
	  Zrmp  (j, i) = sumz;
	}
    }
}

// #######################################################################################################
// Function to calculate angular momentum flux associated with pairs of fully-reconnected solution vectors
// #######################################################################################################
void TJ::GetTorqueFull ()
{
  complex<double> I = complex<double> (0., 1.);

  for (int k = 0; k < nres; k++)
    for (int kp = 0; kp < nres; kp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> Sum = complex<double> (0., 0.);
	  double          q   = Getq (Rgrid[i]);
	  
	  for (int j = 0; j < J; j++)
	    {
	  double mj  = mpol[j];
	  double mnq = mj - ntor * q;
	  
	  Sum += (+ conj(Zf  (j, k, i) + I * Zf  (j, kp, i))
		  * (Psif(j, k, i) + I * Psif(j, kp, i))
		  - conj(Psif(j, k, i) + I * Psif(j, kp, i))
		  * (Zf  (j, k, i) + I * Zf  (j, kp, i))) /mnq;
	    }
	  Sum *= I * M_PI*M_PI * ntor;
	  
	  double torque = real (Sum);
	  
	  Tfull (k, kp, i) = torque;
	}
}

// ###################################################################################################
// Function to calculate angular momentum flux associated with pairs of unreconnected solution vectors
// ###################################################################################################
void TJ::GetTorqueUnrc ()
{
  complex<double> I = complex<double> (0., 1.);

  for (int k = 0; k < nres; k++)
    for (int kp = 0; kp < nres; kp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> Sum = complex<double> (0., 0.);
	  double          q   = Getq (Rgrid[i]);
	  
	  for (int j = 0; j < J; j++)
	    {
	  double mj  = mpol[j];
	  double mnq = mj - ntor * q;
	  
	  Sum += (+ conj(Zu  (j, k, i) + I * Zu  (j, kp, i))
		  * (Psiu(j, k, i) + I * Psiu(j, kp, i))
		  - conj(Psiu(j, k, i) + I * Psiu(j, kp, i))
		  * (Zu  (j, k, i) + I * Zu  (j, kp, i))) /mnq;
	    }
	  Sum *= I * M_PI*M_PI * ntor;
	  
	  double torque = real (Sum);
	  
	  Tunrc (k, kp, i) = torque;
	}
}

// ###############################################
// Function to calculate ideal stabiltiy of plasma
// ###############################################
void TJ::CalculateIdealStability ()
{
  // ---------------
  // Allocate memory
  // ---------------
  Psii.resize(J, J, NDIAG);
  Zi  .resize(J, J, NDIAG);
  Ji  .resize(J, J);
  Wmat.resize(J, J);
  Wher.resize(J, J);
  Want.resize(J, J);
  Wvec.resize(J, J);
  
  Psie.resize(J, J, NDIAG);
  Ze  .resize(J, J, NDIAG);
  Je  .resize(J, J);
  Pres.resize(J, J);

  Wval   = new double[J];
  deltaW = new double[J];
  
  // -----------------------------------------------------
  // Calculate ideal solutions launched from magnetic axis
  // -----------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = YYY(j,   jp, i);
	  complex<double> sumz = YYY(J+j, jp, i);
	  
	  for (int k = 0; k < nres; k++)
	    {
	      sump -= Psiu(j, k, i) * Pia(k, jp);
	      sumz -= Zu  (j, k, i) * Pia(k, jp);
	    }
	  
	  Psii(j, jp, i) = sump;
	  Zi  (j, jp, i) = sumz;
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
	  sum -= Hmat(j, jpp) * Psii(jpp, jp, NDIAG-1);

	Ji(j, jp) = sum;
      }

  // -----------------------------
  // Calculate ideal energy matrix
  // -----------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += conj (Psii(jpp, j, NDIAG-1)) * Ji(jpp, jp);
	
	Wmat(j, jp) = sum;	      
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	Wher (j, jp) = 0.5 * (Wmat (j, jp) + conj (Wmat (jp, j)));
	Want (j, jp) = 0.5 * (Wmat (j, jp) - conj (Wmat (jp, j)));
      }

  // ............................
  // Calculate W-matrix residuals
  // ............................
  double Whmax = 0., Wamax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double whval = abs (Wher (j, jp));
	double waval = abs (Want (j, jp));

	if (whval > Whmax)
	  Whmax = whval;
	if (waval > Wamax)
	  Wamax = waval;	
      }

  printf ("Energy matrix residual: %10.4e\n", Wamax /Whmax);

  // --------------------------------------------------
  // Calculate eigenvalues and eigenvectors of W-matrix
  // --------------------------------------------------
  GetEigenvalues (Wher, Wval, Wvec);

  // -----------------------------
  // Calculate ideal eigenfuctions
  // -----------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      for (int i = 0; i < NDIAG; i++)
	{
	  complex<double> sump = complex<double> (0., 0.);
	  complex<double> sumz = complex<double> (0., 0.);
	  
	  for (int jpp = 0; jpp < J; jpp++)
	    {
	      sump += Psii(j, jpp, i) * Wvec(jpp, jp);
	      sumz += Zi  (j, jpp, i) * Wvec(jpp, jp);
	    }
	  
	  Psie(j, jp, i) = sump;
	  Ze  (j, jp, i) = sumz;
	}

  // ------------------------------
  // Normalize ideal eigenfunctions
  // ------------------------------
  for (int jp = 0; jp < J; jp++)
    {
      double norm = 0.;

      for (int j = 0; j < J; j++)
	norm += real (conj (Psie(j, jp, NDIAG-1)) * Psie(j, jp, NDIAG-1));

      for (int i = 0; i < NDIAG; i++)
	for (int j = 0; j < J; j++)
	  {
	    Psie(j, jp, i) /= sqrt(norm);
	    Ze  (j, jp, i) /= sqrt(norm);
	  }
    }

  // ----------------------------------------------------------------
  // Calculate boundary currents associated with ideal eigenfunctions
  // ----------------------------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = Ze(j, jp, NDIAG-1) /(mpol[j] - ntor*qa);

	for (int jpp = 0; jpp < J; jpp++)
	  sum -= Hmat(j, jpp) * Psie(jpp, jp, NDIAG-1);

	Je(j, jp) = sum;
      }

  // ------------------------
  // Calculate delta-W values
  // ------------------------
  for (int j = 0; j < J; j++)
    {
      double sum = 0.;

      for (int jp = 0; jp < J; jp++)
	sum += real (conj (Psie(jp, j, NDIAG-1)) * Je(jp, j));

      deltaW[j] = sum;
    }

  printf ("Ideal eigenvalues:\n");
  for (int j = 0; j < J; j++)
    printf ("j = %3d  Wval = %10.3e  deltaW = %10.3e\n", j, Wval[j], deltaW[j]);
  
  // --------------------------------------------
  // Check orthonormality of ideal eigenfunctions
  // --------------------------------------------
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> sum = complex<double> (0., 0.);

	for (int jpp = 0; jpp < J; jpp++)
	  sum += conj (Psie(jpp, j, NDIAG-1)) * Je(jpp, jp) /deltaW[jp];

	if (j == jp)
	  Pres(j, jp) = sum - complex<double> (1., 0.);
	else
	  Pres(j, jp) = sum;
      }

  // ---------------------------------
  // Calculate orthogonality residuals
  // ---------------------------------
  double Pmax = 0.;
  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < J; jp++)
      {
	double pval = abs (Pres (j, jp));

	if (pval > Pmax)
	  Pmax = pval;
      }

  printf ("Ideal eigenfunction othogonality residual: %10.4e\n", Pmax);
}

// ##############################################################################
// Function to output visualization data for unreconnected tearing eigenfunctions
// ##############################################################################
void TJ::VisualizeEigenfunctions ()
{
  // ...............
  // Allocate memory
  // ...............
  Psiuf.resize(J,    nres, Nf);
  Zuf  .resize(J,    nres, Nf);
  Psiuv.resize(nres, Nf,   Nw+1);
  Zuv  .resize(nres, Nf,   Nw+1);
  Psirf.resize(J,    Nf);
  Zrf  .resize(J,    Nf);
  Psirv.resize(Nf,   Nw+1);
  Zrv  .resize(Nf,   Nw+1);

  // .........................................................................................
  // Interpolate unreconnected eigenfunction data from diagnostic to visualization radial grid
  // .........................................................................................
  for (int j = 0; j < J; j++)
    for (int k = 0; k < nres; k++)
      {
	// Get data from diagnostic grid
	double* psi_r = new double[NDIAG];
	double* psi_i = new double[NDIAG];
	double* z_r   = new double[NDIAG];
	double* z_i   = new double[NDIAG];

	for (int i = 0; i < NDIAG; i++)
	  {
	    psi_r[i] = real(Psiu(j, k, i));
	    psi_i[i] = imag(Psiu(j, k, i));
	    z_r  [i] = real(Zu  (j, k, i));
	    z_i  [i] = imag(Zu  (j, k, i));
	  }

	// Interpolate data from diagnostic grid
	gsl_interp_accel* psi_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* psi_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* z_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* z_i_acc   = gsl_interp_accel_alloc ();

	gsl_spline* psi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* psi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);

	gsl_spline_init (psi_r_spline, Rgrid, psi_r, NDIAG);
	gsl_spline_init (psi_i_spline, Rgrid, psi_i, NDIAG);
	gsl_spline_init (z_r_spline,   Rgrid, z_r,   NDIAG);
	gsl_spline_init (z_i_spline,   Rgrid, z_i,   NDIAG);

	// Interpolate data onto visualization grid
	for (int i = 0; i < Nf; i++)
	  {
	    double x = gsl_spline_eval (psi_r_spline, rf[i], psi_r_acc);
	    double y = gsl_spline_eval (psi_i_spline, rf[i], psi_i_acc);

	    Psiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (z_r_spline, rf[i], z_r_acc);
	    y = gsl_spline_eval (z_i_spline, rf[i], z_i_acc);

	    Zuf(j, k, i) = complex<double> (x, y);
	  }

	// Clean up
	delete[] psi_r; delete[] psi_i; delete[] z_r; delete[] z_i;

	gsl_spline_free (psi_r_spline);
	gsl_spline_free (psi_i_spline);
	gsl_spline_free (z_r_spline);
	gsl_spline_free (z_i_spline);

	gsl_interp_accel_free (psi_r_acc);
	gsl_interp_accel_free (psi_i_acc);
	gsl_interp_accel_free (z_r_acc);
	gsl_interp_accel_free (z_i_acc);
      }

  // ............................................................
  // Calculate unreconnected eigenfunctions on visualization grid
  // ............................................................
  complex<double> II = complex<double> (0., 1.);

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < Nf; i++)
      for (int l = 0; l <= Nw; l++)
	{
	  double theta = thvals(i, l);
	  double psi_r = 0., psi_i = 0., z_r = 0., z_i = 0.;

	  for (int j = 0; j < J; j++)
	    {
	      double m = mpol[j];

	      psi_r += real (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      psi_i += imag (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      // Omit m=0 component of Z
	      if (MPOL[j] != 0)
		{
		  z_r += real (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		  z_i += imag (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		}
	    }

	  Psiuv(k, i, l) = complex<double> (psi_r, psi_i);
	  Zuv  (k, i, l) = complex<double> (z_r,   z_i);
	}

  // ................................................................................
  // Interpolate ideal RMP response data from diagnostic to visualization radial grid
  // ................................................................................
  for (int j = 0; j < J; j++)
    {
      // Get data from diagnostic grid
      double* psi_r = new double[NDIAG];
      double* psi_i = new double[NDIAG];
      double* z_r   = new double[NDIAG];
      double* z_i   = new double[NDIAG];
      
      for (int i = 0; i < NDIAG; i++)
	{
	  psi_r[i] = real(Psirmp(j, i));
	  psi_i[i] = imag(Psirmp(j, i));
	  z_r  [i] = real(Zrmp  (j, i));
	  z_i  [i] = imag(Zrmp  (j, i));
	}
      
      // Interpolate data from diagnostic grid
      gsl_interp_accel* psi_r_acc = gsl_interp_accel_alloc ();
      gsl_interp_accel* psi_i_acc = gsl_interp_accel_alloc ();
      gsl_interp_accel* z_r_acc   = gsl_interp_accel_alloc ();
      gsl_interp_accel* z_i_acc   = gsl_interp_accel_alloc ();
      
      gsl_spline* psi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
      gsl_spline* psi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
      gsl_spline* z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
      gsl_spline* z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
      
      gsl_spline_init (psi_r_spline, Rgrid, psi_r, NDIAG);
      gsl_spline_init (psi_i_spline, Rgrid, psi_i, NDIAG);
      gsl_spline_init (z_r_spline,   Rgrid, z_r,   NDIAG);
      gsl_spline_init (z_i_spline,   Rgrid, z_i,   NDIAG);
      
      // Interpolate data onto visualization grid
      for (int i = 0; i < Nf; i++)
	{
	  double x = gsl_spline_eval (psi_r_spline, rf[i], psi_r_acc);
	  double y = gsl_spline_eval (psi_i_spline, rf[i], psi_i_acc);
	  
	  Psirf(j, i) = complex<double> (x, y);
	  
	  x = gsl_spline_eval (z_r_spline, rf[i], z_r_acc);
	  y = gsl_spline_eval (z_i_spline, rf[i], z_i_acc);
	  
	  Zrf(j, i) = complex<double> (x, y);
	}
      
      // Clean up
      delete[] psi_r; delete[] psi_i; delete[] z_r; delete[] z_i;
      
      gsl_spline_free (psi_r_spline);
      gsl_spline_free (psi_i_spline);
      gsl_spline_free (z_r_spline);
      gsl_spline_free (z_i_spline);
      
      gsl_interp_accel_free (psi_r_acc);
      gsl_interp_accel_free (psi_i_acc);
      gsl_interp_accel_free (z_r_acc);
      gsl_interp_accel_free (z_i_acc);
    }
  
  // ................................................................
  // Calculate ideal RMP response eigenfunction on visualization grid
  // ................................................................
  for (int i = 0; i < Nf; i++)
    for (int l = 0; l <= Nw; l++)
      {
	double theta = thvals(i, l);
	double psi_r = 0., psi_i = 0., z_r = 0., z_i = 0.;
	
	for (int j = 0; j < J; j++)
	  {
	    double m = mpol[j];
	    
	    psi_r += real (Psirf(j, i) * (cos(m*theta) + II*sin(m*theta)));
	    psi_i += imag (Psirf(j, i) * (cos(m*theta) + II*sin(m*theta)));
	    
	    z_r += real (Zrf(j, i) * (cos(m*theta) + II*sin(m*theta)));
	    z_i += imag (Zrf(j, i) * (cos(m*theta) + II*sin(m*theta)));
	  }
	
	Psirv(i, l) = complex<double> (psi_r, psi_i);
	Zrv  (i, l) = complex<double> (z_r,   z_i);
      }
}

