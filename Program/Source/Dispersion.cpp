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
  Emat .resize(nres, nres);
  Psif .resize(J,    nres, NDIAG);
  Zf   .resize(J,    nres, NDIAG);
  Tf   .resize(nres, NDIAG);
  Psiu .resize(J,    nres, NDIAG);
  Zu   .resize(J,    nres, NDIAG);
  Tu   .resize(nres, NDIAG);
  Tfull.resize(nres, nres, NDIAG);
  Tunrc.resize(nres, nres, NDIAG);

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
	  Ymat(j, jp) = - Zs(j, jp) /(mpol[j] - ntor*qa);
	else
	  Ymat(j, jp) = complex<double> (0., 0.);

	for (int k = 0; k < J; k++)
	  Ymat(j, jp) += Hmat(j, k) * Psis(k, jp);
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
	Fmat(j, jp) = Pis(j, jp);

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
  // Calculate fully-reconnected tearing eigenfunctions
  // ..................................................
  for (int i = 0; i < NDIAG; i++)
    {
      for (int j = 0; j < J; j++)
	for (int k = 0; k < nres; k++)
	  {
	    Psif(j, k, i) = YYY(j, J+k, i);

	    for (int jp = 0; jp < J; jp++)
	      Psif(j, k, i) += YYY(j, jp, i) * Omat(jp, k);

	    Zf(j, k, i) = YYY(J+j, J+k, i);

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

  // ........................................................
  // Calculate unreconnected eigenfunction visualization data
  // ........................................................
  VisualizeEigenfunctions ();
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

  // ..........................................................................................
  // Interpolate unreconnected eigenfunction data from diagnostic to visulalization radial grid
  // ..........................................................................................
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
		  z_r += real (Zuf  (j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		  z_i += imag (Zuf  (j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		}
	    }

	  Psiuv(k, i, l) = complex<double> (psi_r, psi_i);
	  Zuv  (k, i, l) = complex<double> (z_r,   z_i);
	}
}

