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
  psiu .resize(J,    nres, NDIAG);
  zu   .resize(J,    nres, NDIAG);
  chiu .resize(J,    nres, NDIAG);
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
	if (FREE < 0)
	  Xmat(j, jp) = complex<double> (0., 0.);
	else
	  Xmat(j, jp) = Za(j, jp) /(mpol[j] - ntor*qa);

	if (FREE > 0 || FREE < 0)
	  for (int k = 0; k < J; k++)
	    Xmat(j, jp) -= Hmat(j, k) * Psia(k, jp);
	else
	   for (int k = 0; k < J; k++)
	    Xmat(j, jp) -= Gmat(j, k) * Psia(k, jp);
      }

  for (int j = 0; j < J; j++)
    for (int jp = 0; jp < nres; jp++)
      {
	if (FREE < 0)
	  Ymat(j, jp) = complex<double> (0., 0.);
	else
	  Ymat(j, jp) = - Zs(j, jp) /(mpol[j] - ntor*qa) /dPi(jp);

	if (FREE > 0 || FREE < 0)
	  for (int k = 0; k < J; k++)
	    Ymat(j, jp) += Hmat(j, k) * Psis(k, jp) /dPi(jp);
	else
	   for (int k = 0; k < J; k++)
	     Ymat(j, jp) += Gmat(j, k) * Psis(k, jp) /dPi(jp);
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

  if (FVAL)
    {
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
  printf ("E-matrix Hermitian test residual: %10.3e\n", Eerr);

  // ..............................................
  // Calculate unreconnected tearing eigenfunctions
  // ..............................................
  for (int i = 0; i < NDIAG; i++)
    {
      double r = Rgrid[i];
      double q = Getq (r);
	    
      for (int j = 0; j < J; j++)
	for (int k = 0; k < nres; k++)
	  {
	    Psiu(j, k, i) = complex<double> (0., 0.);

	    for (int kp = 0; kp < nres; kp++)
	      Psiu(j, k, i) += Psif(j, kp, i) * Emat(kp, k);

	    Zu(j, k, i) = complex<double> (0., 0.);

	    for (int kp = 0; kp < nres; kp++)
	      Zu(j, k, i) += Zf(j, kp, i) * Emat(kp, k);

	    double delta = double (mres[k]) * sres[k] * ISLAND /2. /rres[k];
	    double mnq   = double (mres[k]) - ntor * q;
	    double imnq  = mnq / (delta*delta + mnq*mnq);

	    double PSI;
	    if (mres[k] == 1)
	      PSI = epsa * (rres[k] * gres[k] * sres[k] /qres[k]) * ISLAND /abs (Emat(k, k));
	    else
	      PSI = epsa * (ISLAND*ISLAND /16.) *  (gres[k] * sres[k] /qres[k]); 

	    psiu(j, k, i) = PSI * Psiu(j, k, i); 
	    
	    if (MPOL[j] == mres[k])
	      zu(j, k, i) = PSI * (Zu(j, k, i) + Getkm (r, MPOL[j]) * Psiu(j, k, i)) * imnq;
	    else
	      zu(j, k, i) = PSI * (Zu(j, k, i) + Getkm (r, MPOL[j]) * Psiu(j, k, i)) /(mpol[j] - ntor * q);
	  }
    }

  for (int i = 0; i < NDIAG; i++)
    {
      double r = Rgrid[i];
      double q = Getq (r);
      
      Array<complex<double>,2> LLmmp (J, J);
      Array<complex<double>,2> MMmmp (J, J);
      Array<complex<double>,2> NNmmp (J, J);
      Array<complex<double>,2> PPmmp (J, J);
      
      GetMatrices (r, LLmmp, MMmmp, NNmmp, PPmmp);

      for (int j = 0; j < J; j++)
	for (int k = 0; k < nres; k++)
	  {
	    double delta = double (mres[k]) * sres[k] * ISLAND /2. /rres[k];
	    double mnq   = double (mres[k]) - ntor * q;
	    double imnq  = mnq / (delta*delta + mnq*mnq);

	    double PSI;
	    if (mres[k] == 1)
	      PSI = epsa * (rres[k] * gres[k] * sres[k] /qres[k]) * ISLAND /abs (Emat(k, k));
	    else
	      PSI = epsa * (ISLAND*ISLAND /16.) *  (gres[k] * sres[k] /qres[k]); 
	    
	    if (MPOL[j] == 0)
	      chiu(j, k, i) = complex<double> (0., 0.);
	    else
	      {
		chiu(j, k, i) = - epsa*epsa * ntor*ntor * r*r * zu(j, k, i);

		for (int jp = 0; jp < J; jp++)
		  {
		    if (MPOL[jp] == mres[k])
		      chiu(j, k, i) += (LLmmp(j, jp) * Zu(jp, k, i) + MMmmp(j, jp) * Psiu(jp, k, i)) * imnq;
		    else
		      chiu(j, k, i) += (LLmmp(j, jp) * Zu(jp, k, i) + MMmmp(j, jp) * Psiu(jp, k, i)) /(mpol[jp] - ntor * q);
		  }

		chiu(j, k, i) *= PSI /mpol[j];
	      }
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

