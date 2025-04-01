// Visulalize.cpp

#include "TJ.h"

// ##############################################################################
// Function to output visualization data for unreconnected tearing eigenfunctions
// ##############################################################################
void TJ::VisualizeEigenfunctions ()
{
  // ...............
  // Allocate memory
  // ...............
  Psiuf.resize(J, nres, Nf);
  Zuf  .resize(J, nres, Nf);
  psiuf.resize(J, nres, Nf);
  zuf  .resize(J, nres, Nf);
  chiuf.resize(J, nres, Nf);
  xiuf .resize(J, nres, Nf);

  Psiuv.resize(nres, Nf, Nw+1);
  Zuv  .resize(nres, Nf, Nw+1);
  zuv  .resize(nres, Nf, Nw+1);
  chiuv.resize(nres, Nf, Nw+1);
  bRc  .resize(nres, Nf, Nw+1);
  bRs  .resize(nres, Nf, Nw+1);
  bZc  .resize(nres, Nf, Nw+1);
  bZs  .resize(nres, Nf, Nw+1);
  bPc  .resize(nres, Nf, Nw+1);
  bPs  .resize(nres, Nf, Nw+1);
  xic  .resize(nres, Nf, Nw+1);
  xis  .resize(nres, Nf, Nw+1);
 
  // .........................................................................................
  // Interpolate unreconnected eigenfunction data from diagnostic to visualization radial grid
  // .........................................................................................
  for (int j = 0; j < J; j++)
    for (int k = 0; k < nres; k++)
      {
	// Get data from diagnostic grid
	double* psi_r = new double[NDIAG];
	double* psi_i = new double[NDIAG];
	double* Z_r   = new double[NDIAG];
	double* Z_i   = new double[NDIAG];

	double* p_r   = new double[NDIAG];
	double* p_i   = new double[NDIAG];
	double* z_r   = new double[NDIAG];
	double* z_i   = new double[NDIAG];
	double* chi_r = new double[NDIAG];
	double* chi_i = new double[NDIAG];
	double* xi_r  = new double[NDIAG];
	double* xi_i  = new double[NDIAG];

	for (int i = 0; i < NDIAG; i++)
	  {
	    psi_r[i] = real (Psiu(j, k, i));
	    psi_i[i] = imag (Psiu(j, k, i));
	    Z_r  [i] = real (Zu  (j, k, i));
	    Z_i  [i] = imag (Zu  (j, k, i));

	    p_r  [i] = real (psiu(j, k, i));
	    p_i  [i] = imag (psiu(j, k, i));
	    z_r  [i] = real (zu  (j, k, i));
	    z_i  [i] = imag (zu  (j, k, i));
	    chi_r[i] = real (chiu(j, k, i));
	    chi_i[i] = imag (chiu(j, k, i));
	    xi_r [i] = real (xiu (j, k, i));
	    xi_i [i] = imag (xiu (j, k, i));
	  }

	// Interpolate data from diagnostic grid
	gsl_interp_accel* psi_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* psi_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* Z_i_acc   = gsl_interp_accel_alloc ();
	
	gsl_interp_accel* p_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* p_i_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* z_r_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* z_i_acc   = gsl_interp_accel_alloc ();
	gsl_interp_accel* chi_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* chi_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* xi_r_acc  = gsl_interp_accel_alloc ();
	gsl_interp_accel* xi_i_acc  = gsl_interp_accel_alloc ();
	
	gsl_spline* psi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* psi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	
	gsl_spline* p_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* p_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* z_r_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* z_i_spline   = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* chi_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* chi_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* xi_r_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* xi_i_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);

	gsl_spline_init (psi_r_spline, Rgrid, psi_r, NDIAG);
	gsl_spline_init (psi_i_spline, Rgrid, psi_i, NDIAG);
	gsl_spline_init (Z_r_spline,   Rgrid, Z_r,   NDIAG);
	gsl_spline_init (Z_i_spline,   Rgrid, Z_i,   NDIAG);

	gsl_spline_init (p_r_spline,   Rgrid, p_r,   NDIAG);
	gsl_spline_init (p_i_spline,   Rgrid, p_i,   NDIAG);
	gsl_spline_init (z_r_spline,   Rgrid, z_r,   NDIAG);
	gsl_spline_init (z_i_spline,   Rgrid, z_i,   NDIAG);
	gsl_spline_init (chi_r_spline, Rgrid, chi_r, NDIAG);
	gsl_spline_init (chi_i_spline, Rgrid, chi_i, NDIAG);
	gsl_spline_init (xi_r_spline,  Rgrid, xi_r,  NDIAG);
	gsl_spline_init (xi_i_spline,  Rgrid, xi_i,  NDIAG);
		
	// Interpolate data onto visualization grid
	for (int i = 0; i < Nf; i++)
	  {
	    double x = gsl_spline_eval (psi_r_spline, rf[i], psi_r_acc);
	    double y = gsl_spline_eval (psi_i_spline, rf[i], psi_i_acc);

	    Psiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (Z_r_spline, rf[i], Z_r_acc);
	    y = gsl_spline_eval (Z_i_spline, rf[i], Z_i_acc);

	    Zuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (p_r_spline, rf[i], p_r_acc);
	    y = gsl_spline_eval (p_i_spline, rf[i], p_i_acc);

	    psiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (z_r_spline, rf[i], z_r_acc);
	    y = gsl_spline_eval (z_i_spline, rf[i], z_i_acc);

	    zuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (chi_r_spline, rf[i], chi_r_acc);
	    y = gsl_spline_eval (chi_i_spline, rf[i], chi_i_acc);

	    chiuf(j, k, i) = complex<double> (x, y);

	    x = gsl_spline_eval (xi_r_spline, rf[i], xi_r_acc);
	    y = gsl_spline_eval (xi_i_spline, rf[i], xi_i_acc);

	    xiuf(j, k, i) = complex<double> (x, y);
	  }

	// Clean up
	delete[] psi_r; delete[] psi_i; delete[] Z_r;  delete[] Z_i;
	delete[] chi_r; delete[] chi_i; delete[] z_r;  delete[] z_i;
	delete[] p_r;   delete[] p_i;   delete[] xi_r; delete[] xi_i;

	gsl_spline_free (psi_r_spline);
	gsl_spline_free (psi_i_spline);
	gsl_spline_free (Z_r_spline);
	gsl_spline_free (Z_i_spline);
	gsl_spline_free (p_r_spline);
	gsl_spline_free (p_i_spline);
	gsl_spline_free (z_r_spline);
	gsl_spline_free (z_i_spline);
	gsl_spline_free (chi_r_spline);
	gsl_spline_free (chi_i_spline);
	gsl_spline_free (xi_r_spline);
	gsl_spline_free (xi_i_spline);
	
	gsl_interp_accel_free (psi_r_acc);
	gsl_interp_accel_free (psi_i_acc);
	gsl_interp_accel_free (Z_r_acc);
	gsl_interp_accel_free (Z_i_acc);
	gsl_interp_accel_free (p_r_acc);
	gsl_interp_accel_free (p_i_acc);
	gsl_interp_accel_free (z_r_acc);
	gsl_interp_accel_free (z_i_acc);
	gsl_interp_accel_free (chi_r_acc);
	gsl_interp_accel_free (chi_i_acc);
	gsl_interp_accel_free (xi_r_acc);
	gsl_interp_accel_free (xi_i_acc);
      }

  // ............................................................
  // Calculate unreconnected eigenfunctions on visualization grid
  // ............................................................
  complex<double> II = complex<double> (0., 1.);

  for (int k = 0; k < nres; k++)
    for (int i = 0; i < Nf; i++)
      for (int l = 0; l <= Nw; l++)
	{
	  double R     = RR    (i, l);
	  double r     = rvals (i, l);
	  double theta = thvals(i, l);
	  double psi_r = 0., psi_i = 0., Z_r   = 0., Z_i   = 0.;
	  double z_r   = 0., z_i   = 0., chi_r = 0., chi_i = 0.;
	  double bR_c  = 0., bR_s  = 0., bZ_c  = 0., bZ_s  = 0., bP_c = 0., bP_s = 0.;
	  double xi_c  = 0., xi_s  = 0.;

	  for (int j = 0; j < J; j++)
	    {
	      double m = mpol[j];

	      psi_r += real (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      psi_i += imag (Psiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      // Omit m=0 component of Z
	      if (MPOL[j] != 0)
		{
		  Z_r += real (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		  Z_i += imag (Zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
		}

	      z_r += real (zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      z_i += imag (zuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      chi_r += real (chiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));
	      chi_i += imag (chiuf(j, k, i) * (cos(m*theta) + II*sin(m*theta)));

	      bR_c +=   dRdr(i, l) * (   real (psiuf(j, k, i)) * sin(m*theta) + imag (psiuf(j, k, i)) * cos(m*theta))
	 	      + dRdt(i, l) * (   real (chiuf(j, k, i)) * cos(m*theta) - imag (chiuf(j, k, i)) * sin(m*theta));
	      bR_s +=   dRdr(i, l) * ( - real (psiuf(j, k, i)) * cos(m*theta) + imag (psiuf(j, k, i)) * sin(m*theta))
		      + dRdt(i, l) * (   real (chiuf(j, k, i)) * sin(m*theta) + imag (chiuf(j, k, i)) * cos(m*theta));
	      bZ_c +=   dZdr(i, l) * (   real (psiuf(j, k, i)) * sin(m*theta) + imag (psiuf(j, k, i)) * cos(m*theta))
	 	      + dZdt(i, l) * (   real (chiuf(j, k, i)) * cos(m*theta) - imag (chiuf(j, k, i)) * sin(m*theta));
	      bZ_s +=   dZdr(i, l) * ( - real (psiuf(j, k, i)) * cos(m*theta) + imag (psiuf(j, k, i)) * sin(m*theta))
		      + dZdt(i, l) * (   real (chiuf(j, k, i)) * sin(m*theta) + imag (chiuf(j, k, i)) * cos(m*theta));
	      bP_c +=                    real (zuf  (j, k, i)) * cos(m*theta) - imag (zuf  (j, k, i)) * sin(m*theta);
	      bP_s +=                    real (zuf  (j, k, i)) * sin(m*theta) + imag (zuf  (j, k, i)) * cos(m*theta);

	      xi_c +=                    real (xiuf (j, k, i)) * cos(m*theta) - imag (xiuf (j, k, i)) * sin(m*theta);
	      xi_s +=                    real (xiuf (j, k, i)) * sin(m*theta) + imag (xiuf (j, k, i)) * cos(m*theta);
	    }

	  Psiuv(k, i, l) = complex<double> (psi_r, psi_i);
	  Zuv  (k, i, l) = complex<double> (Z_r,   Z_i);
	  zuv  (k, i, l) = complex<double> (z_r,   z_i);
	  chiuv(k, i, l) = complex<double> (chi_r, chi_i);

	  bRc(k, i, l) = - B0 * bR_c /r/R/R;
	  bRs(k, i, l) = - B0 * bR_s /r/R/R;
	  bZc(k, i, l) = - B0 * bZ_c /r/R/R;
	  bZs(k, i, l) = - B0 * bZ_s /r/R/R;
	  bPc(k, i, l) =   B0 * epsa * ntor * bP_c /R;
	  bPs(k, i, l) =   B0 * epsa * ntor * bP_s /R;

	  xic(k, i, l) = xi_c;
	  xis(k, i, l) = xi_s;
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
