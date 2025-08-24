// Visualize.cpp

#include "Vertical.h"

// ##############################################################################
// Function to output visualization data for unreconnected tearing eigenfunctions
// ##############################################################################
void Vertical::VisualizeEigenfunctions ()
{
  // ...............
  // Allocate memory
  // ...............
  Psiuf.resize(J, K, Nf);
  Zuf  .resize(J, K, Nf);
  Psiuv.resize(K, Nf, Nw+1);
  Zuv  .resize(K, Nf, Nw+1);
 
  // .........................................................................................
  // Interpolate no-wall ideal eigenfunction data from diagnostic to visualization radial grid
  // .........................................................................................
  for (int j = 0; j < J; j++)
    for (int k = 0; k < K; k++)
      {
	// Get data from diagnostic grid
	double* psi_r = new double[NDIAG];
	double* psi_i = new double[NDIAG];
	double* Z_r   = new double[NDIAG];
	double* Z_i   = new double[NDIAG];
		
	for (int i = 0; i < NDIAG; i++)
	  {
	    psi_r[i] = real (Psie(j, k, i));
	    psi_i[i] = imag (Psie(j, k, i));
	    Z_r  [i] = real (Ze  (j, k, i));
	    Z_i  [i] = imag (Ze  (j, k, i));
	  }

	// Interpolate data from diagnostic grid
	gsl_interp_accel* psi_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* psi_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_i_acc   = gsl_interp_accel_alloc ();
	
	gsl_spline* psi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* psi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	
	gsl_spline_init (psi_r_spline, Rgrid, psi_r, NDIAG);
	gsl_spline_init (psi_i_spline, Rgrid, psi_i, NDIAG);
	gsl_spline_init (Z_r_spline,   Rgrid, Z_r,   NDIAG);
	gsl_spline_init (Z_i_spline,   Rgrid, Z_i,   NDIAG);

	// Interpolate data onto visualization grid
	for (int i = 0; i < Nf; i++)
	  {
	    double x = gsl_spline_eval (psi_r_spline, rf[i], psi_r_acc);
	    double y = gsl_spline_eval (psi_i_spline, rf[i], psi_i_acc);

	    Psiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (Z_r_spline, rf[i], Z_r_acc);
	    y = gsl_spline_eval (Z_i_spline, rf[i], Z_i_acc);

	    Zuf(j, k, i) = complex<double> (x, y);
	  }
	
	// Clean up
	delete[] psi_r; delete[] psi_i; delete[] Z_r;   delete[] Z_i;

	gsl_spline_free (psi_r_spline);
	gsl_spline_free (psi_i_spline);
	gsl_spline_free (Z_r_spline);
	gsl_spline_free (Z_i_spline);
      	
	gsl_interp_accel_free (psi_r_acc);
	gsl_interp_accel_free (psi_i_acc);
	gsl_interp_accel_free (Z_r_acc);
	gsl_interp_accel_free (Z_i_acc);
      }
	
  // ......................................................
  // Calculate no-wall eigenfunctions on visualization grid
  // ......................................................
  complex<double> II = complex<double> (0., 1.);

  for (int k = 0; k < K; k++)
    for (int i = 0; i < Nf; i++)
      for (int l = 0; l <= Nw; l++)
	{
	  double R     = RR    (i, l);
	  double r     = rvals (i, l);
	  double theta = thvals(i, l);
	  double psi_r = 0., psi_i = 0., Z_r = 0., Z_i = 0.;

	  for (int j = 0; j < J; j++)
	    {
	      double m = mpol[j];

	      psi_r += real (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      psi_i += imag (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      Z_r += real (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      Z_i += imag (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	    }
	   
	  Psiuv(k, i, l) = complex<double> (psi_r, psi_i);
	  Zuv  (k, i, l) = complex<double> (Z_r,   Z_i);
	}

  // ...............
  // Allocate memory
  // ...............
  pPsiuf.resize(J, K, Nf);
  pZuf  .resize(J, K, Nf);
  pPsiuv.resize(K, Nf, Nw+1);
  pZuv  .resize(K, Nf, Nw+1);
 
  // ..............................................................................................
  // Interpolate perfect-wall ideal eigenfunction data from diagnostic to visualization radial grid
  // ..............................................................................................
  for (int j = 0; j < J; j++)
    for (int k = 0; k < K; k++)
      {
	// Get data from diagnostic grid
	double* psi_r = new double[NDIAG];
	double* psi_i = new double[NDIAG];
	double* Z_r   = new double[NDIAG];
	double* Z_i   = new double[NDIAG];
		
	for (int i = 0; i < NDIAG; i++)
	  {
	    psi_r[i] = real (pPsie(j, k, i));
	    psi_i[i] = imag (pPsie(j, k, i));
	    Z_r  [i] = real (pZe  (j, k, i));
	    Z_i  [i] = imag (pZe  (j, k, i));
	  }

	// Interpolate data from diagnostic grid
	gsl_interp_accel* psi_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* psi_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_i_acc   = gsl_interp_accel_alloc ();
	
	gsl_spline* psi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* psi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	
	gsl_spline_init (psi_r_spline, Rgrid, psi_r, NDIAG);
	gsl_spline_init (psi_i_spline, Rgrid, psi_i, NDIAG);
	gsl_spline_init (Z_r_spline,   Rgrid, Z_r,   NDIAG);
	gsl_spline_init (Z_i_spline,   Rgrid, Z_i,   NDIAG);

	// Interpolate data onto visualization grid
	for (int i = 0; i < Nf; i++)
	  {
	    double x = gsl_spline_eval (psi_r_spline, rf[i], psi_r_acc);
	    double y = gsl_spline_eval (psi_i_spline, rf[i], psi_i_acc);

	    pPsiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (Z_r_spline, rf[i], Z_r_acc);
	    y = gsl_spline_eval (Z_i_spline, rf[i], Z_i_acc);

	    pZuf(j, k, i) = complex<double> (x, y);
	  }
	
	// Clean up
	delete[] psi_r; delete[] psi_i; delete[] Z_r;   delete[] Z_i;

	gsl_spline_free (psi_r_spline);
	gsl_spline_free (psi_i_spline);
	gsl_spline_free (Z_r_spline);
	gsl_spline_free (Z_i_spline);
      	
	gsl_interp_accel_free (psi_r_acc);
	gsl_interp_accel_free (psi_i_acc);
	gsl_interp_accel_free (Z_r_acc);
	gsl_interp_accel_free (Z_i_acc);
      }
	
  // ...........................................................
  // Calculate perfect-wall eigenfunctions on visualization grid
  // ...........................................................
  for (int k = 0; k < K; k++)
    for (int i = 0; i < Nf; i++)
      for (int l = 0; l <= Nw; l++)
	{
	  double R     = RR    (i, l);
	  double r     = rvals (i, l);
	  double theta = thvals(i, l);
	  double psi_r = 0., psi_i = 0., Z_r = 0., Z_i = 0.;

	  for (int j = 0; j < J; j++)
	    {
	      double m = mpol[j];

	      psi_r += real (pPsiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      psi_i += imag (pPsiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      Z_r += real (pZuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      Z_i += imag (pZuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	    }
	   
	  pPsiuv(k, i, l) = complex<double> (psi_r, psi_i);
	  pZuv  (k, i, l) = complex<double> (Z_r,   Z_i);
	}
 }

