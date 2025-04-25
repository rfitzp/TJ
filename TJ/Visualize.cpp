// Visualize.cpp

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

  neuf .resize(J+1, nres, Nf);
  Teuf .resize(J+1, nres, Nf);
  dneuf.resize(J+1, nres, Nf);
  dTeuf.resize(J+1, nres, Nf);

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

  nec .resize(nres, Nf, Nw+1, NPHI);
  Tec .resize(nres, Nf, Nw+1, NPHI);
  dnec.resize(nres, Nf, Nw+1, NPHI);
  dTec.resize(nres, Nf, Nw+1, NPHI);

  bReqc .resize(nres, 2*Nf, NPHI);
  neeqc .resize(nres, 2*Nf, NPHI);
  Teeqc .resize(nres, 2*Nf, NPHI);
  dneeqc.resize(nres, 2*Nf, NPHI);
  dTeeqc.resize(nres, 2*Nf, NPHI);
  Teeqd.resize (nres, 2*Nf, NPHI);
  dTeeqd.resize(nres, 2*Nf, NPHI);

  Leq  = new double[2*Nf]; 
  Lres = new double[nres];
  Rres = new double[nres];
  Ores = new double[nres];
  Xres = new double[nres];
  PP   = new double[NPHI];

  for (int i = 0; i < NPHI; i++)
    PP[i] = double (i) * 2.*M_PI /double (NPHI);

  itheta       = new double           [2*Nf];
  Teeqcspline  = new gsl_spline*      [nres*NPHI];
  dTeeqcspline = new gsl_spline*      [nres*NPHI];
  Teeqcacc     = new gsl_interp_accel*[nres*NPHI];
  dTeeqcacc    = new gsl_interp_accel*[nres*NPHI];

  for (int k = 0; k < nres; k++)
    for (int np = 0; np < NPHI; np++)
      {
	Teeqcspline [k*NPHI + np] = gsl_spline_alloc (gsl_interp_cspline, 2*Nf);
	dTeeqcspline[k*NPHI + np] = gsl_spline_alloc (gsl_interp_cspline, 2*Nf);

	Teeqcacc    [k*NPHI + np] = gsl_interp_accel_alloc ();
	dTeeqcacc   [k*NPHI + np] = gsl_interp_accel_alloc ();
      }
 
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

	double* ne_r  = new double[NDIAG];
	double* ne_i  = new double[NDIAG];
	double* Te_r  = new double[NDIAG];
	double* Te_i  = new double[NDIAG];
	double* dne_r = new double[NDIAG];
	double* dne_i = new double[NDIAG];
	double* dTe_r = new double[NDIAG];
	double* dTe_i = new double[NDIAG];
		
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

	    if (TEMP)
	      {
		ne_r [i] = real (neu (j, k, i));
		ne_i [i] = imag (neu (j, k, i));
		Te_r [i] = real (Teu (j, k, i));
		Te_i [i] = imag (Teu (j, k, i));
		dne_r[i] = real (dneu(j, k, i));
		dne_i[i] = imag (dneu(j, k, i));
		dTe_r[i] = real (dTeu(j, k, i));
		dTe_i[i] = imag (dTeu(j, k, i));
	      }
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

	gsl_interp_accel* ne_r_acc  = gsl_interp_accel_alloc ();
	gsl_interp_accel* ne_i_acc  = gsl_interp_accel_alloc ();
	gsl_interp_accel* Te_r_acc  = gsl_interp_accel_alloc ();
	gsl_interp_accel* Te_i_acc  = gsl_interp_accel_alloc ();
	gsl_interp_accel* dne_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* dne_i_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* dTe_r_acc = gsl_interp_accel_alloc ();
	gsl_interp_accel* dTe_i_acc = gsl_interp_accel_alloc ();
      
	
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

	gsl_spline* ne_r_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* ne_i_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Te_r_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* Te_i_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* dne_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* dne_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* dTe_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	gsl_spline* dTe_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
      	
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

	if (TEMP)
	  {
	    gsl_spline_init (ne_r_spline,  Rgrid, ne_r,  NDIAG);
	    gsl_spline_init (ne_i_spline,  Rgrid, ne_i,  NDIAG);
	    gsl_spline_init (Te_r_spline,  Rgrid, Te_r,  NDIAG);
	    gsl_spline_init (Te_i_spline,  Rgrid, Te_i,  NDIAG);		
	    gsl_spline_init (dne_r_spline, Rgrid, dne_r, NDIAG);
	    gsl_spline_init (dne_i_spline, Rgrid, dne_i, NDIAG);
	    gsl_spline_init (dTe_r_spline, Rgrid, dTe_r, NDIAG);
	    gsl_spline_init (dTe_i_spline, Rgrid, dTe_i, NDIAG);
	  }
		
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

	    if (TEMP)
	      {
		x = gsl_spline_eval (ne_r_spline, rf[i], ne_r_acc);
		y = gsl_spline_eval (ne_i_spline, rf[i], ne_i_acc);

		neuf(j, k, i) = complex<double> (x, y);
		
		x = gsl_spline_eval (Te_r_spline, rf[i], Te_r_acc);
		y = gsl_spline_eval (Te_i_spline, rf[i], Te_i_acc);
		
		Teuf(j, k, i) = complex<double> (x, y);
		
		x = gsl_spline_eval (dne_r_spline, rf[i], dne_r_acc);
		y = gsl_spline_eval (dne_i_spline, rf[i], dne_i_acc);
		
		dneuf(j, k, i) = complex<double> (x, y);
		
		x = gsl_spline_eval (dTe_r_spline, rf[i], dTe_r_acc);
		y = gsl_spline_eval (dTe_i_spline, rf[i], dTe_i_acc);
		
		dTeuf(j, k, i) = complex<double> (x, y);
	      }
	  }

	// Clean up
	delete[] psi_r; delete[] psi_i; delete[] Z_r;   delete[] Z_i;
	delete[] chi_r; delete[] chi_i; delete[] z_r;   delete[] z_i;
	delete[] p_r;   delete[] p_i;   delete[] xi_r;  delete[] xi_i;

	delete[] dne_r; delete[] dne_i; delete[] dTe_r; delete[] dTe_i; 
	delete[] ne_r;  delete[] ne_i;  delete[] Te_r;  delete[] Te_i;

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

	gsl_spline_free (ne_r_spline);
	gsl_spline_free (ne_i_spline);
	gsl_spline_free (Te_r_spline);
	gsl_spline_free (Te_i_spline);
	gsl_spline_free (dne_r_spline);
	gsl_spline_free (dne_i_spline);
	gsl_spline_free (dTe_r_spline);
	gsl_spline_free (dTe_i_spline);
      	
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

	gsl_interp_accel_free (ne_r_acc);
	gsl_interp_accel_free (ne_i_acc);
	gsl_interp_accel_free (Te_r_acc);
	gsl_interp_accel_free (Te_i_acc);
	gsl_interp_accel_free (dne_r_acc);
	gsl_interp_accel_free (dne_i_acc);
	gsl_interp_accel_free (dTe_r_acc);
	gsl_interp_accel_free (dTe_i_acc);
      }

  if (TEMP)
    {
      for (int k = 0; k < nres; k++)
	{
	  // Get data from diagnostic grid
	  double* ne_r  = new double[NDIAG];
	  double* ne_i  = new double[NDIAG];
	  double* Te_r  = new double[NDIAG];
	  double* Te_i  = new double[NDIAG];
	  double* dne_r = new double[NDIAG];
	  double* dne_i = new double[NDIAG];
	  double* dTe_r = new double[NDIAG];
	  double* dTe_i = new double[NDIAG];
	  
	  for (int i = 0; i < NDIAG; i++)
	    {
	      ne_r [i] = real (neu (J, k, i));
	      ne_i [i] = imag (neu (J, k, i));
	      Te_r [i] = real (Teu (J, k, i));
	      Te_i [i] = imag (Teu (J, k, i));
	      dne_r[i] = real (dneu(J, k, i));
	      dne_i[i] = imag (dneu(J, k, i));
	      dTe_r[i] = real (dTeu(J, k, i));
	      dTe_i[i] = imag (dTeu(J, k, i));
	    }
	  
	  // Interpolate data from diagnostic grid
	  gsl_interp_accel* ne_r_acc  = gsl_interp_accel_alloc ();
	  gsl_interp_accel* ne_i_acc  = gsl_interp_accel_alloc ();
	  gsl_interp_accel* Te_r_acc  = gsl_interp_accel_alloc ();
	  gsl_interp_accel* Te_i_acc  = gsl_interp_accel_alloc ();
	  gsl_interp_accel* dne_r_acc = gsl_interp_accel_alloc ();
	  gsl_interp_accel* dne_i_acc = gsl_interp_accel_alloc ();
	  gsl_interp_accel* dTe_r_acc = gsl_interp_accel_alloc ();
	  gsl_interp_accel* dTe_i_acc = gsl_interp_accel_alloc ();
	  
	  gsl_spline* ne_r_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* ne_i_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* Te_r_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* Te_i_spline  = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* dne_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* dne_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* dTe_r_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  gsl_spline* dTe_i_spline = gsl_spline_alloc (gsl_interp_cspline, NDIAG);
	  
	  gsl_spline_init (ne_r_spline,  Rgrid, ne_r,  NDIAG);
	  gsl_spline_init (ne_i_spline,  Rgrid, ne_i,  NDIAG);
	  gsl_spline_init (Te_r_spline,  Rgrid, Te_r,  NDIAG);
	  gsl_spline_init (Te_i_spline,  Rgrid, Te_i,  NDIAG);		
	  gsl_spline_init (dne_r_spline, Rgrid, dne_r, NDIAG);
	  gsl_spline_init (dne_i_spline, Rgrid, dne_i, NDIAG);
	  gsl_spline_init (dTe_r_spline, Rgrid, dTe_r, NDIAG);
	  gsl_spline_init (dTe_i_spline, Rgrid, dTe_i, NDIAG);
	  
	  // Interpolate data onto visualization grid
	  for (int i = 0; i < Nf; i++)
	    {
	      double x = gsl_spline_eval (ne_r_spline, rf[i], ne_r_acc);
	      double y = gsl_spline_eval (ne_i_spline, rf[i], ne_i_acc);
	      
	      neuf(J, k, i) = complex<double> (x, y);
	      
	      x = gsl_spline_eval (Te_r_spline, rf[i], Te_r_acc);
	      y = gsl_spline_eval (Te_i_spline, rf[i], Te_i_acc);
	      
	      Teuf(J, k, i) = complex<double> (x, y);
	      
	      x = gsl_spline_eval (dne_r_spline, rf[i], dne_r_acc);
	      y = gsl_spline_eval (dne_i_spline, rf[i], dne_i_acc);
	      
	      dneuf(J, k, i) = complex<double> (x, y);
	      
	      x = gsl_spline_eval (dTe_r_spline, rf[i], dTe_r_acc);
	      y = gsl_spline_eval (dTe_i_spline, rf[i], dTe_i_acc);
	      
	      dTeuf(J, k, i) = complex<double> (x, y);
	    }
	  
	  // Clean up
	  delete[] dne_r; delete[] dne_i; delete[] dTe_r; delete[] dTe_i; 
	  delete[] ne_r;  delete[] ne_i;  delete[] Te_r;  delete[] Te_i; 
	  
	  gsl_spline_free (ne_r_spline);
	  gsl_spline_free (ne_i_spline);
	  gsl_spline_free (Te_r_spline);
	  gsl_spline_free (Te_i_spline);
	  gsl_spline_free (dne_r_spline);
	  gsl_spline_free (dne_i_spline);
	  gsl_spline_free (dTe_r_spline);
	  gsl_spline_free (dTe_i_spline);
      
	  gsl_interp_accel_free (ne_r_acc);
	  gsl_interp_accel_free (ne_i_acc);
	  gsl_interp_accel_free (Te_r_acc);
	  gsl_interp_accel_free (Te_i_acc);
	  gsl_interp_accel_free (dne_r_acc);
	  gsl_interp_accel_free (dne_i_acc);
	  gsl_interp_accel_free (dTe_r_acc);
	  gsl_interp_accel_free (dTe_i_acc);
	}
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
	  double bR_c  = 0., bR_s  = 0., bZ_c  = 0., bZ_s  = 0.;
	  double bP_c  = 0., bP_s  = 0., xi_c  = 0., xi_s  = 0.;

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

	      xi_c  +=                   real (xiuf (j, k, i)) * cos(m*theta) - imag (xiuf (j, k, i)) * sin(m*theta);
	      xi_s  +=                   real (xiuf (j, k, i)) * sin(m*theta) + imag (xiuf (j, k, i)) * cos(m*theta);
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

	  xic (k, i, l) = xi_c;
	  xis (k, i, l) = xi_s;
	}

  if (TEMP)
    {
      for (int k = 0; k < nres; k++)
	for (int i = 0; i < Nf; i++)
	  for (int l = 0; l <= Nw; l++)
	    for (int np = 0; np < NPHI; np++)
	      {
		double island;
		if (k > ISLAND.size () - 1)
		  island = 0.;
		else
		  island = ISLAND[k];
		
		double r     = rvals (i, l);
		double theta = thvals(i, l);
		
		double ne_c  = real (neuf(J, k, i));
		double Te_c  = real (Teuf(J, k, i));		       
		double dne_c = 0.;
		double dTe_c = 0.;
		
		for (int j = 0; j < J; j++)
		  {
		    double m = mpol[j];
		    
		    complex<double> eik = cos (m*theta - ntor*PP[np]) + II * sin (m*theta - ntor*PP[np]);
		    
		    ne_c  += real (neuf (j, k, i) * eik);
		    Te_c  += real (Teuf (j, k, i) * eik);
		    dne_c += real (dneuf(j, k, i) * eik);
		    dTe_c += real (dTeuf(j, k, i) * eik);
		    
		    double x = r - rres[k];
		    
		    if (fabs(x) < island && m == mres[k])
		      {
			for (int nh = 2; nh < Nh; nh++)
			  {
			    complex<double> eikh = cos (double (nh) * (m*theta - ntor*PP[np])) + II * sin (double (nh) * (m*theta - ntor*PP[np]));
			    
			    double dTh = gsl_spline_eval (dThspline[k*Nh + nh], x /island, dThacc[k*Nh + nh]);
				
			    ne_c  += real (nepres[k] * island * dTh * eikh);
			    Te_c  += real (Tepres[k] * island * dTh * eikh);
			    dne_c += real (nepres[k] * island * dTh * eikh);
			    dTe_c += real (Tepres[k] * island * dTh * eikh);
			  }
		      }
		  }
		
		nec (k, i, l, np) = ne_c;
		Tec (k, i, l, np) = Te_c;
		dnec(k, i, l, np) = dne_c;
		dTec(k, i, l, np) = dTe_c;
	      }
      
      // ............................................
      // Calculate quantities on tilted central chord
      // ............................................
      double wt = tilt * M_PI/180;
      for (int i = 0; i < Nf; i++)
	Leq[i] = (Req[i] - 1.) /cos(wt);

      for (int k = 0; k < nres; k++)
	for (int i = 0; i < Nf; i++)
	  for (int np = 0; np < NPHI; np++)
	    {
	      double island;
	      if (k > ISLAND.size () - 1)
		island = 0.;
	      else
		island = ISLAND[k];
		
	      double R     = Req   [i];
	      double r     = req   [i];
	      double theta = teq   [i];
	      double dRdr  = dRdreq[i];
	      double dRdt  = dRdteq[i];
	      double dZdr  = dZdreq[i];
	      double dZdt  = dZdteq[i];
	      double bR_c  = 0., bZ_c  = 0., dne_c = 0., dTe_c = 0.;
	      
	      double ne_c = real (neuf(J, k, Nf-1-i));
	      double Te_c = real (Teuf(J, k, Nf-1-i));
	      
	      for (int j = 0; j < J; j++)
		{
		  double m = mpol[j];
		  
		  complex<double> eik = cos (m*theta - ntor*PP[np]) + II * sin (m*theta - ntor*PP[np]);
		  
		  bR_c += real ((dRdr * II * psiuf(j, k, Nf-1-i) - dRdt * chiuf(j, k, Nf-1-i)) * eik);
		  bZ_c += real ((dZdr * II * psiuf(j, k, Nf-1-i) - dZdt * chiuf(j, k, Nf-1-i)) * eik);
		  
		  ne_c  += real (neuf (j, k, Nf-1-i) * eik);
		  Te_c  += real (Teuf (j, k, Nf-1-i) * eik);
		  dne_c += real (dneuf(j, k, Nf-1-i) * eik);
		  dTe_c += real (dTeuf(j, k, Nf-1-i) * eik);
		  
		  double x = r - rres[k];
		  
		  if (fabs(x) < island && m == mres[k])
		    {
		      for (int nh = 2; nh < Nh; nh++)
			{
			  complex<double> eikh = cos (double (nh) * (m*theta - ntor*PP[np])) + II * sin (double (nh) * (m*theta - ntor*PP[np]));
			  
			  double dTh = gsl_spline_eval (dThspline[k*Nh + nh], x /island, dThacc[k*Nh + nh]);
			  
			  ne_c  += real (nepres[k] * island * dTh * eikh);
			  Te_c  += real (Tepres[k] * island * dTh * eikh);
			  dne_c += real (nepres[k] * island * dTh * eikh);
			  dTe_c += real (Tepres[k] * island * dTh * eikh);
			  
			}
		    }
		}
	      
	      bReqc (k, i, np) = - B0 * (bR_c * cos(wt) - bZ_c * sin(wt)) /r/R/R;
	      neeqc (k, i, np) = ne_c;
	      Teeqc (k, i, np) = Te_c;
	      dneeqc(k, i, np) = dne_c;
	      dTeeqc(k, i, np) = dTe_c;
	    }
      
      wt = tilt * M_PI/180;
      for (int i = 0; i < Nf; i++)
	Leq[Nf+i] = (Req[Nf+i] - 1.) /cos(wt);
      
      double* data1 = new double[Nf];
      double* data2 = new double[Nf];

      for (int i = 0; i < Nf; i++)
	{
	  data1[i] = req[Nf+i];
	  data2[i] = Leq[Nf+i];
	}

      gsl_spline*       dataspline = gsl_spline_alloc (gsl_interp_cspline, Nf);
      gsl_interp_accel* dataacc    = gsl_interp_accel_alloc ();

      gsl_spline_init (dataspline, data1, data2, Nf);

      for (int k = 0; k < nres; k++)
	{
	  Lres[k] = gsl_spline_eval (dataspline, rres[k], dataacc);
	}

      for (int i = 0; i < Nf; i++)
	{
	  data1[i] = req[Nf+i];
	  data2[i] = Req[Nf+i];
	}

      gsl_spline_init (dataspline, data1, data2, Nf);

      for (int k = 0; k < nres; k++)
	{
	  Rres[k] = gsl_spline_eval (dataspline, rres[k],                               dataacc);
	  Ores[k] = gsl_spline_eval (dataspline, rres[k] - delta[k]*ISLAND[k]/sqrt(8.), dataacc);
	  Xres[k] = gsl_spline_eval (dataspline, rres[k] + delta[k]*ISLAND[k]/sqrt(8.), dataacc);
	}

      delete[] data1; delete[] data2;
      gsl_spline_free (dataspline);
      gsl_interp_accel_free (dataacc);
      
      for (int k = 0; k < nres; k++)
	for (int i = 0; i < Nf; i++)
	  for (int np = 0; np < NPHI; np++)
	    {
	      double island;
	      if (k > ISLAND.size () - 1)
		island = 0.;
	      else
		island = ISLAND[k];
	      
	      double R     = Req   [Nf+i];
	      double r     = req   [Nf+i];
	      double theta = teq   [Nf+i];
	      double dRdr  = dRdreq[Nf+i];
	      double dRdt  = dRdteq[Nf+i];
	      double dZdr  = dZdreq[Nf+i];
	      double dZdt  = dZdteq[Nf+i];
	      double bR_c  = 0., bZ_c  = 0., dne_c = 0., dTe_c = 0.;
	      
	      double ne_c = real (neuf(J, k, i));
	      double Te_c = real (Teuf(J, k, i));
	      
	      for (int j = 0; j < J; j++)
		{
		  double m = mpol[j];
		  
		  complex<double> eik = cos (m*theta - ntor*PP[np]) + II * sin (m*theta - ntor*PP[np]);
		  
		  bR_c += real ((dRdr * II * psiuf(j, k, i) - dRdt * chiuf(j, k, i)) * eik);
		  bZ_c += real ((dZdr * II * psiuf(j, k, i) - dZdt * chiuf(j, k, i)) * eik);
		  
		  ne_c  += real (neuf (j, k, i) * eik);
		  Te_c  += real (Teuf (j, k, i) * eik);
		  dne_c += real (dneuf(j, k, i) * eik);
		  dTe_c += real (dTeuf(j, k, i) * eik);
		  
		  double x = r - rres[k];
		  
		  if (fabs(x) < island && m == mres[k])
		    {
		      for (int nh = 2; nh < Nh; nh++)
			{
			  complex<double> eikh = cos (double (nh) * (m*theta - ntor*PP[np])) + II * sin (double (nh) * (m*theta - ntor*PP[np]));
			  
			  double dTh = gsl_spline_eval (dThspline[k*Nh + nh], x /island, dThacc[k*Nh + nh]);
			  
			  ne_c  += real (nepres[k] * island * dTh * eikh);
			  Te_c  += real (Tepres[k] * island * dTh * eikh);
			  dne_c += real (nepres[k] * island * dTh * eikh);
			  dTe_c += real (Tepres[k] * island * dTh * eikh);
			  
			}
		    }
		}
	      
	      bReqc (k, Nf+i, np) = - B0 * (bR_c * cos(wt) - bZ_c * sin(wt)) /r/R/R;
	      neeqc (k, Nf+i, np) = ne_c;
	      Teeqc (k, Nf+i, np) = Te_c;
	      dneeqc(k, Nf+i, np) = dne_c;
	      dTeeqc(k, Nf+i, np) = dTe_c;
	    }

      // ###########################################################
      // Correct temperature signals for relativistic ece broadening
      // ###########################################################
      double m_e = 9.1093837015e-31;
      double c   = 2.997924580e8;
      double e   = 1.602176634e-19;

      for (int i = 0; i < 2*Nf; i++)
	itheta[i] = m_e*c*c /e /Teeq[i];

      for (int k = 0; k < nres; k++)
	for (int np = 0; np < NPHI; np++)
	  {
	    double* data1 = new double[2*Nf];
	    double* data2 = new double[2*Nf];

	    for (int i = 0; i < 2*Nf; i++)
	      {
		data1[i] = Teeqc (k, i, np);
		data2[i] = dTeeqc(k, i, np);
	      }

	    gsl_spline_init (Teeqcspline [k*NPHI + np], Req, data1, 2*Nf);
	    gsl_spline_init (dTeeqcspline[k*NPHI + np], Req, data2, 2*Nf);

	    delete[] data1; delete[] data2; 
	  }

      double h, t_err;
      int    rept, count;

      double  R;
      double* Y    = new double[2*nres*NPHI];
      double* err  = new double[2*nres*NPHI];

      for (int i = 1; i < 2*Nf; i++)
	{
	  iomega = i;

	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		dTeeqd(k, 0, np);
	      }
	  
	  rhs_chooser = 0;
	  count       = 0;
	  h           = h0;
	  R           = Req[0];
	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		Y[            k*NPHI + np] = 0.;
		Y[nres*NPHI + k*NPHI + np] = 0.;
	      }

	  do
	    {
	      CashKarp45Adaptive (2*nres*NPHI, R, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (R < Req[i]);
	  CashKarp45Fixed (2*nres*NPHI, R, Y, err, Req[i] - R);

	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		dTeeqd(k, i, np) = Y[k*NPHI + np] /Y[nres*NPHI + k*NPHI + np];
	      }

	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		Teeqd(k, 0, np);
	      }
	  
	  rhs_chooser = 1;
	  count       = 0;
	  h           = h0;
	  R           = Req[0];
	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		Y[            k*NPHI + np] = 0.;
		Y[nres*NPHI + k*NPHI + np] = 0.;
	      }

	  do
	    {
	      CashKarp45Adaptive (2*nres*NPHI, R, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (R < Req[i]);
	  CashKarp45Fixed (2*nres*NPHI, R, Y, err, Req[i] - R);

	  for (int k = 0; k < nres; k++)
	    for (int np = 0; np < NPHI; np++)
	      {
		Teeqd(k, i, np) = Y[k*NPHI + np] /Y[nres*NPHI + k*NPHI + np];
	      }
	}

      delete[] Y; delete[] err;
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

// #####################################
// Evaluate right-hand sides of ece odes
// #####################################
void TJ::CashKarp45Rhs (double R, double* Y, double* dYdR)
{
  double Romega = Req   [iomega];
  double ith    = itheta[iomega];
  
  for (int k = 0; k < nres; k++)
    for (int np = 0; np < NPHI; np++)
      {
	double dTe;
	if (R > Req[2*Nf-1])
	  {
	    dTe = dTeeqc(k, 2*Nf-1, np);
	  }
	else
	  {
	    if (rhs_chooser == 0)
	      dTe = gsl_spline_eval (dTeeqcspline[k*NPHI + np], R, dTeeqcacc[k*NPHI + np]);
	    else
	      dTe = gsl_spline_eval (Teeqcspline [k*NPHI + np], R, Teeqcacc [k*NPHI + np]);
	  }

	dYdR[            k*NPHI + np] =  dTe * (1. - R*R /Romega*Romega) * exp (- ith * (Romega /R - 1.));
	dYdR[nres*NPHI + k*NPHI + np] =        (1. - R*R /Romega*Romega) * exp (- ith * (Romega /R - 1.));
      }
}
