// RMP.cpp

#include "TJ.h"

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

  Psirmps = new complex<double>[Nw+1];
  Psixs   = new complex<double>[Nw+1];

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

  // ..............................................
  // Calculate Psi_x and Psi_rmp on plasma boundary
  // ..............................................
  for (int i = 0; i <= Nw; i++)
    {
      double theta = tbound[i];
      
      complex<double> sump = complex<double> (0., 0.);
      complex<double> sumx = complex<double> (0., 0.);
      
      for (int j = 0; j < J; j++)
	{
	  sump += Psirmp(j, NDIAG-1) * complex<double> (cos (mpol[j] * theta), sin (mpol[j] * theta));
	  sumx += Psix[j]            * complex<double> (cos (mpol[j] * theta), sin (mpol[j] * theta));
	}
      
      Psirmps[i] = sump;
      Psixs  [i] = sumx;
    }
 }

// #################################################################################
// Function to output visualization data for resonant magnetic perturbation response
// #################################################################################
void TJ::VisualizeRMP ()
{
  // ...............
  // Allocate memory
  // ...............
  Psirf.resize(J,  Nf);
  Zrf  .resize(J,  Nf);
  Psirv.resize(Nf, Nw+1);
  Zrv  .resize(Nf, Nw+1);
  
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
  complex<double> II = complex<double> (0., 1.);

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
