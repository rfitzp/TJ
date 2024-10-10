// ODESolve.cpp

#include "TJ.h"

// ###################################
// Function to solve outer region odes
// ###################################
void TJ::ODESolve ()
{
  printf ("ODE solution data:\n");
  
  // +++++++++++++++
  // Allocate memory
  // +++++++++++++++
  K = J + nres;
  Array<complex<double>,2> YY (2*J, K);
  YYY.resize (2*J, K, NDIAG);
  for (int i = 0; i < 2*J; i++)
    for (int j = 0; j < K; j++)
      for (int k = 0; k < NDIAG; k++)
	YYY(i, j, k) = complex<double> (0., 0.);

  Pi.resize (nres, K);
  for (int i = 0; i < nres; i++)
    for (int j = 0; j < K; j++)
      Pi(i, j) = complex<double> (0., 0.);

  dPi.resize (nres);
  for (int i = 0; i < nres; i++)
       dPi(i) = complex<double> (0., 0.);

  Rgrid = new double[NDIAG];
  for (int i = 0; i < NDIAG; i++)
    Rgrid[i] = EPS + (1. - EPS) * double (i) /double (NDIAG - 1);

  Pgrid = new double[NDIAG];
  for (int i = 0; i < NDIAG; i++)
    Pgrid[i] = gsl_spline_eval (Pspline, Rgrid[i], Pacc);

  hode = new double[NDIAG];
  eode = new double[NDIAG];

  Ttest.resize (K, NDIAG);
  Pnorm.resize (K, NDIAG);
  Znorm.resize (K, NDIAG);

  for (int i = 0; i < NDIAG; i++)
    {
      hode[i] = h0;
      eode[i] = acc;
      for (int j = 0; j < K; j++)
	{
	  Ttest (j, i) = 0.;
	  Pnorm (j, i) = 0.;
	  Znorm (j, i) = 0.;
	}
    }

  // ++++++++++++++++++++++++++++++++++++++++++++
  // Initialize solution vectors at magnetic axis
  // ++++++++++++++++++++++++++++++++++++++++++++
  double rr = EPS;
  LaunchAxis (rr, YY);

  printf ("Initializing solution vectors at magnetic axis:             r = %11.4e: Torque test = %11.4e\n", rr, TorqueTest (rr, YY));

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Integrate solution vectors launched from magnetic axis to first rational surface
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  double rx = rres[0] - DEL;
  int    nf = (NFIX == 0) ? 0 : int ((rx) * double (NFIX));
  SegmentFixup (rr, rx, nf, YY);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++
  // Jump solution vectors across first rational surface
  // +++++++++++++++++++++++++++++++++++++++++++++++++++
  printf ("Integrating solution vectors to rational surface (%3d/%3d): r = %11.4e: Torque test = %11.4e\n", 1, nres, rr, TorqueTest (rr, YY));
  JumpRational   (0, rr, YY);
  LaunchRational (0, rr, YY);
  printf ("Integrating solution vectors to rational surface (%3d/%3d): r = %11.4e: Torque test = %11.4e\n", 1, nres, rr, TorqueTest (rr, YY));

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Integrate solution vectors across other rational surfaces
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int j = 1; j < nres; j++)
    {
      rx = rres[j] - DEL;
      nf = (NFIX == 0) ? 0 : int ((rx - rres[j-1]) * double (NFIX));
      SegmentFixup (rr, rx, nf, YY);
      
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Jump solution vectors across other rational surfaces
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++
      printf ("Integrating solution vectors to rational surface (%3d/%3d): r = %11.4e: Torque test = %11.4e\n", j+1, nres, rr, TorqueTest (rr, YY));
      JumpRational   (j, rr, YY);
      LaunchRational (j, rr, YY);
      printf ("Integrating solution vectors to rational surface (%3d/%3d): r = %11.4e: Torque test = %11.4e\n", j+1, nres, rr, TorqueTest (rr, YY));
    }

  // +++++++++++++++++++++++++++++++++++++++++++++
  // Integrate solution vectors to plasma boundary
  // +++++++++++++++++++++++++++++++++++++++++++++
  rx = 1.;
  nf = (NFIX == 0) ? 0 : int ((rx - rres[nres-1]) * double (NFIX));
  SegmentFixup (rr, rx, nf, YY);
  
  printf ("Integrating solution vectors to plasma boundary:            r = %11.4e: Torque test = %11.4e\n", rx, TorqueTest (rx, YY));

  int index = NDIAG - 1;

  for (int j = 0; j < 2*J; j++)
    for (int jp = 0; jp < K; jp++)
      YYY(j, jp, index) = YY(j, jp);

  hode[index] = hode[index-1];
  eode[index] = eode[index-1];
  for (int i = 0; i < K; i++)
    {
      double Pnm, Znm;
      GetNorms (YY, i, Pnm, Znm);
      
      Pnorm(i, index) = Pnm;
      Znorm(i, index) = Znm;
    }

  // ..................................................
  // Calculate torques associated with solution vectors
  // ..................................................
  for (int i = 0; i < NDIAG; i++)
    for (int jp = 0; jp < K; jp++)
      {
	complex<double> I   = complex<double> (0., 1.);
	complex<double> Sum = complex<double> (0., 0.);
	double          q   = Getq (Rgrid[i]);
	
	for (int j = 0; j < J; j++)
	  {
	    double mj  = mpol[j];
	    double mnq = mj - ntor * q;
	    
	    Sum += (conj(YYY(J+j, jp, i)) * YYY(j, jp, i) - conj(YYY(j, jp, i)) * YYY(J+j, jp, i)) /mnq;
	  }
	Sum *= I * M_PI*M_PI * ntor;

	Ttest(jp, i) = real(Sum);
      }
}

// ######################################################
// Function to launch solution vectors from magnetic axis
// ######################################################
void TJ::LaunchAxis (double r, Array<complex<double>,2> YY)
{
  Array<complex<double>,2> PPsi (J, K);
  Array<complex<double>,2> ZZ   (J, K);

  double ppp = Getppp (r);
  double q   = Getq   (r);
  double s2  = Gets2  (r);

  double* Hnp = new double[Ns+1];
  double* Vnp = new double[Ns+1];

  for (int n = 1; n <= Ns; n++)
    Hnp[n] = GetHnp (n, r);
  for (int n = 2; n <= Ns; n++)
    Vnp[n] = GetVnp (n, r);
   
  // Launch J independent solutions vectors from magnetic axis
  for (int j = 0; j < J; j++)
    {
      int    MJ  = MPOL[j];
      double mj  = mpol[j];
      double mnq = mj - ntor * q;

      for (int k = 0; k < J; k++)
	{
	  PPsi(k, j) = complex<double> (0., 0.);
	  ZZ  (k, j) = complex<double> (0., 0.);
	}

       if (MJ == 0)
	 ZZ(j, j) = complex<double> (1., 0.);
       else if (MJ > 0)
	 {
	   PPsi(j, j) = complex<double> (1., 0.);
	   ZZ  (j, j) = + mnq /mj;
	 }
       else if (MJ < 0)
	 {
	   PPsi(j, j) = complex<double> (1., 0.);
	   ZZ  (j, j) = - mnq /mj;
	 }
    }

  // Launch nres null solution vectors from magnetic axis
  for (int j = 0; j < nres; j++)
    {
      for (int k = 0; k < J; k++)
	{
	  PPsi(k, J+j) = complex<double> (0., 0.);
	  ZZ  (k, J+j) = complex<double> (0., 0.);
	}
    }
  
  PackYY (PPsi, ZZ, YY);

  int index = 0;
 
  for (int j = 0; j < 2*J; j++)
    for (int jp = 0; jp < K; jp++)
      YYY(j, jp, index) = YY(j, jp);
  
  hode[index] = h0;
  eode[index] = acc;
 
  for (int i = 0; i < K; i++)
    {
      double Pnm, Znm;
      GetNorms (YY, i, Pnm, Znm);

      Pnorm(i, index) = Pnm;
      Znorm(i, index) = Znm;
    }

  delete[] Hnp; delete[] Vnp;
}

// ###############################################################
// Function to launch solution vector from jresth rational surface
// ###############################################################
void TJ::LaunchRational (int jres, double r, Array<complex<double>,2> YY)
{
  // ...............................
  // Set rational surface quantities
  // ...............................
  double mm = double(mres[jres]);
  double rm = rres[jres];
  double qm = qres[jres];
  int    km = Jres[jres];

  // ........................
  // Calculate L1, P1, and T1
  // ........................
  double sm, L1, P1, T1;
  GetL1P1T1 (rm, mm, km, sm, L1, P1, T1);

  // ..........................
  // Calculate Mercier indicies
  // ..........................
  Array<complex<double>,2> Lmat(J, J), Mmat(J, J), Nmat(J, J), Pmat(J, J);
  GetMatrices (rm, Lmat, Mmat, Nmat, Pmat);

  double L0 = - real (Lmat(km, km)) /mm /sm;
  double P0 = - real (Pmat(km, km)) /mm /sm;
  double RT =   sqrt (0.25 + L0 * P0);

  double nuL = 0.5 - RT;
  double nuS = 0.5 + RT;

  double bS   = nuS /L0;
  double dnuS = pow (DEL, nuS);
  double norm = pow (rm, nuS) * sqrt ((nuS - nuL) /real(Lmat(km, km)));

  // ..................
  // Calculate akt, bkt
  // ..................
  complex<double>* akt = new complex<double>[J];
  complex<double>* bkt = new complex<double>[J];
  for (int k = 0; k < J; k++)
    {
      if (k != km)
	{
	  akt[k] = - (Lmat(k, km) /L0  + Mmat(k, km) /nuS) /mm /sm;
	  bkt[k] = - (Pmat(k, km) /nuS + Nmat(k, km) /L0 ) /mm /sm;
	}
      else
	{
	  akt[k] = complex<double> (0., 0.);
	  bkt[k] = complex<double> (0., 0.);
	}
    }

  Array<complex<double>,2> PPsi(J, K), ZZ(J, K);

  UnpackYY (YY, PPsi, ZZ);

  for (int k = 0; k < J; k++)
    {
      if (k != km)
	{
	  PPsi(k, J+jres) += akt[k]*dnuS /norm;
	  ZZ  (k, J+jres) += bkt[k]*dnuS /norm;
	}
      else
	{
	  PPsi(k, J+jres) += dnuS      /norm;
	  ZZ  (k, J+jres) += dnuS * bS /norm;
	}
    }

  PackYY (PPsi, ZZ, YY);
  dPi (jres) = complex<double> (1., 0.);

  delete[] akt; delete[] bkt;
}

// ##########################################################################################################
// Function to integrate multiple solution vectors from given value of r to r = rx while performing nf fixups
// ##########################################################################################################
void TJ::SegmentFixup (double& r, double rx, int nf, Array<complex<double>,2> YY)
{
  double rw   = r, RE;
  int    lmax = (nf == 0) ? 1 : nf;

  for (int l = 1; l <= lmax; l++)
    {
      // Determine next fixup radius
      if (nf < 2)
	RE = rx;
      else
	RE = exp(log(rw) + (log(rx) - log(rw)) * double(l) /double(nf));
      
      // Integrate solution vectors to next fixup radius
      printf ("Integrating solution vectors to fixup radius (%3d/%3d):     r = %11.4e", l, nf, RE);
      Segment (r, RE, YY);
      
      // Fixup solution vector
      if (nf > 0)
	{
	  Fixup (r, YY);
	}

      printf (": Torque test = %11.4e\n", TorqueTest (r, YY));
    }
}

// ##################################################################################
// Function to integrate multiple solution vectors from given value of r to to r = rx
// ##################################################################################
void TJ::Segment (double& r, double rx, Array<complex<double>,2> YY)
{
  int              neqns = (2*J)*(K);
  double           h, t_err;
  int              rept;
  complex<double>* Y   = new complex<double>[neqns];
  complex<double>* err = new complex<double>[neqns];

  h     = h0;
  count = 0;
  
  int index = 0;
  for (int i = 0; i < NDIAG; i++)
    if (Rgrid[i] < r)
      index = i;

  if (Rgrid[index+1] < rx)
    {
      do
	{
	  index++; 
	  
	  UnpackYY (YY, Y);
	  do
	    {
	      CashKarp45Adaptive (neqns, r, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (r < Rgrid[index]);
	  CashKarp45Fixed (neqns, r, Y, err, Rgrid[index] - r);
	  PackYY (Y, YY);
	  
	  for (int j = 0; j < 2*J; j++)
	    for (int jp = 0; jp < K; jp++)
	      YYY(j, jp, index) = YY(j, jp);
	  
	  hode[index] = h;
	  eode[index] = t_err;
	  for (int i = 0; i < K; i++)
	    {
	      double Pnm, Znm;
	      GetNorms (YY, i, Pnm, Znm);
	      
	      Pnorm(i, index) = Pnm;
	      Znorm(i, index) = Znm;
	    }
	}
      while (Rgrid[index+1] < rx);
    }
  
  UnpackYY (YY, Y);
  do
    {
      CashKarp45Adaptive (neqns, r, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
   }
  while (r < rx);
  CashKarp45Fixed (neqns, r, Y, err, rx - r);
  PackYY (Y, YY);

  delete[] Y; delete[] err;
}

// ######################################################
// Function to perform fixup of multiple solution vectors
// ######################################################
void TJ::Fixup (double r, Array<complex<double>,2> YY)
{
  // +++++++++++++++++++++++++++++++++++
  // Find index of m=0 poloidal harmonic
  // +++++++++++++++++++++++++++++++++++
  int i0 = -1;
  for (int i = 0; i < J; i++)
    if (MPOL[i] == 0)
      i0 = i;

  // +++++++++++++++++++++++++++++++++++++++++++++++++++
  // Range of poloidal mode numbers does not include m=0
  // (so all mode numbers are positive)
  // +++++++++++++++++++++++++++++++++++++++++++++++++++
  if (i0 == -1)
    {
      for (int i = J-1; i >= 0; i--)
	{
	  complex<double> yii = YY(i, i);
	  
	  for (int j = i-1; j >= 0; j--)
	    {
	      complex<double> yij = YY(i, j);

	      // Fixup solution vectors
	      for (int k = 0; k < 2*J; k++)
		{
		  complex<double> yki = YY(k, i);
		  complex<double> ykj = YY(k, j);

		  YY(k, j) =  ykj - (yij /yii) * yki;
		}

	      // Fixup previously calculated solution vectors
	      for (int l = 0; l < NDIAG; l++)
		{
		  if (r > Rgrid[l])
		    {
		      for (int k = 0; k < 2*J; k++)
			{
			  complex<double> Yki = YYY(k, i, l);
			  complex<double> Ykj = YYY(k, j, l);
			  
			  YYY(k, j, l) = Ykj - (yij /yii) * Yki;
			}
		    }
		}

	      // Fixup Pi matrix
	      for (int k = 0; k < nres; k++)
		{
		  if (r > rres[k])
		    Pi(k, j) = Pi(k, j) - (yij /yii) * Pi(k, i);
		}
	    }
	}
    }
  // +++++++++++++++++++++++++++++++++++++++++++++++
  // Range of poloidal mode numbers does include m=0 
  // +++++++++++++++++++++++++++++++++++++++++++++++
  else
    {
      // Range down positive dominant mode numbers
      for (int i = J-1; i > i0; i--)
	{
	  complex<double> yii = YY(i, i);
	  
	  for (int j = i-1; j >= 0; j--)
	    {
	      complex<double> yij = YY(i, j);

	      // Fixup solution vectors
	      for (int k = 0; k < 2*J; k++)
		{
		  complex<double> yki = YY(k, i);
		  complex<double> ykj = YY(k, j);

		  YY(k, j) = ykj - (yij /yii) * yki;
		}

	      // Fixup previously calculated solution vectors
	      for (int l = 0; l < NDIAG; l++)
		{
		  if (r > Rgrid[l])
		    {
		      for (int k = 0; k < 2*J; k++)
			{
			  complex<double> Yki = YYY(k, i, l);
			  complex<double> Ykj = YYY(k, j, l);
			  
			  YYY(k, j, l) = Ykj - (yij /yii) * Yki;
			}
		    }
		}
	      
	      // Fixup Pi matrix
	      for (int k = 0; k < nres; k++)
		{
		  if (r > rres[k])
		    Pi(k, j) = Pi(k, j) - (yij /yii) * Pi(k, i);
		}
	    }
	}
      // Range up negative dominant mode numbers
      for (int i = 0; i < i0; i++)
	{
	  complex<double> yii = YY(i, i);

	  for (int j = i+1; j < J; j++)
	    {
	      complex<double> yij = YY(i, j);

	      // Fixup solution vectors
	      for (int k = 0; k < 2*J; k++)
		{
		  complex<double> yki = YY(k, i);
		  complex<double> ykj = YY(k, j);

		  YY(k, j) = ykj - (yij /yii) * yki;
		}

	      // Fixup previously calculated solution vectors
	      for (int l = 0; l < NDIAG; l++)
		{
		  if (r > Rgrid[l])
		    {
		      for (int k = 0; k < 2*J; k++)
			{
			  complex<double> Yki = YYY(k, i, l);
			  complex<double> Ykj = YYY(k, j, l);
		      
			  YYY(k, j, l) = Ykj - (yij /yii) * Yki;
			}
		    }
		}
	      
	      // Fixup Pi matrix
	      for (int k = 0; k < nres; k++)
		{
		  if (r > rres[k])
		    Pi(k, j) = Pi(k, j) - (yij /yii) * Pi(k, i);
		}
	    }
	}
    }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Renormalize independent solution vectors launched from magnetic axis
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int i = 0; i < J; i++)
    {
      complex<double> yii = (MPOL[i] == 0) ? YY(J+i, i) : YY(i, i);

      if (abs (yii) > 1.e-15)
	{
          // Renormalize solution vectors
	  for (int k = 0; k < 2*J; k++)
	    {
	      complex<double> yki = YY(k, i);
	      
	      YY(k, i) = yki /yii;
	      
	      for (int l = 0; l < NDIAG; l++)
		{
		  if (r > Rgrid[l])
		    {
		      complex<double> Yki = YYY(k, i, l);
		      
		      YYY(k, i, l) = Yki /yii;
		    }
		}
	    }
	  
	  // Renormalize Pi matrix
	  for (int k = 0; k < nres; k++)
	    {
	      if (r > rres[k])
		{
		  complex<double> pki = Pi(k, i);
		  
		  Pi(k, i) = pki /yii;
		}
	    }
	}
    }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Renormalize independent solution vectors launched from rational surfaces
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int i = 0; i < nres; i++)
    {
      complex<double> yii = YY(Jres[i], J+i);

      if (abs (yii) > 1.e-15)
	{
          // Renormalize solution vectors
	  for (int k = 0; k < 2*J; k++)
	    {
	      complex<double> yki = YY(k, J+i);
	      
	      YY(k, J+i) = yki /yii;
	      
	      for (int l = 0; l < NDIAG; l++)
		{
		  if (r > Rgrid[l])
		    {
		      complex<double> Yki = YYY(k, J+i, l);
		      
		      YYY(k, J+i, l) = Yki /yii;
		    }
		}
	    }
	  
	  // Renormalize Pi matrix
	  for (int k = 0; k < nres; k++)
	    {
	      if (r > rres[k])
		{
		  complex<double> pki = Pi(k, J+i);
		  
		  Pi(k, J+i) = pki /yii;
		}
	    }

	  // Renormalize dPi vector
	  complex<double> dpi = dPi(i);
	  
	  dPi(i) = dpi /yii;
	}
    }
}

// ##########################################################
// Function to evaluate right-hand sides of outer region odes
// ##########################################################
void TJ::Rhs (double r, complex<double>* Y, complex<double>* dYdr)
{
  Array<complex<double>,2> PPsi    (J, K);
  Array<complex<double>,2> ZZ      (J, K);
  Array<complex<double>,2> dPPsidr (J, K);
  Array<complex<double>,2> dZZdr   (J, K);
  Array<complex<double>,2> LLmmp   (J, J);
  Array<complex<double>,2> MMmmp   (J, J);
  Array<complex<double>,2> NNmmp   (J, J);
  Array<complex<double>,2> PPmmp   (J, J);

  GetMatrices (r, LLmmp, MMmmp, NNmmp, PPmmp);
  
  UnpackY (Y, PPsi, ZZ);

  double q   = Getq (r);
  double s   = Gets (r);
  double qp  = q * s /r;
  double nq  = ntor * q;
  double nqp = ntor * qp;

  for (int i = 0; i < K; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  double mjmnq = mpol[j] - nq;
	
	  dPPsidr(j, i) = complex<double> (0., 0.);
	  dZZdr  (j, i) = - nqp * ZZ(j, i) /mjmnq;

	  for (int jp = 0; jp < J; jp++)
	    {
	      double mjpmnq = mpol[jp] - nq;
	      
	      dPPsidr(j, i) += (LLmmp(j, jp) * ZZ(jp, i) + MMmmp(j, jp) * PPsi(jp, i)) /mjpmnq /r;
	      dZZdr  (j, i) += (NNmmp(j, jp) * ZZ(jp, i) + PPmmp(j, jp) * PPsi(jp, i)) /mjpmnq /r;
	    }
	}
    }
  
  PackY (dPPsidr, dZZdr, dYdr);
}

// ###################################
// Function to pack YY solution vector
// ###################################
 void TJ::PackYY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, Array<complex<double>,2> YY)
{
  for (int i = 0; i < K; i++)
    for (int j = 0; j < J; j++)
      {
	YY(j,   i) = PPsi(j, i);
	YY(J+j, i) = ZZ  (j, i); 
      }
}

// ######################################
// Function to unpack YY solution vectors 
// ######################################
void TJ::UnpackYY (Array<complex<double>,2> YY, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ)
{
  for (int i = 0; i < K; i++)
    for (int j = 0; j < J; j++)
      {
	PPsi(j, i) = YY(j,   i);
	ZZ  (j, i) = YY(J+j, i);
      }
}

// ####################################
// Function to pack YY solution vectors 
// ####################################
void TJ::PackYY (complex<double>* Y, Array<complex<double>,2> YY)
{
  int index = 0;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < 2*J; j++)
      {
	YY(j, i) = Y[index]; index++;
      }
}

// #####################################
// Function to unpack YY solution vector
// #####################################
void TJ::UnpackYY (Array<complex<double>,2> YY, complex<double>* Y)
{
  int index = 0;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < 2*J; j++)
      {
	Y[index] = YY(j, i); index++;
      }
}

// ##################################
// Function to pack Y solution vector
// ##################################
void TJ::PackY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, complex<double>* Y)
{
  int index = 0;
  for (int i = 0; i < K; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  Y[index] = PPsi(j, i); index++;
	}
      for (int j = 0; j < J; j++)
	{
	  Y[index] = ZZ(j, i); index++;
	}
    }
}

// ####################################
// Function to unpack Y solution vector
// ####################################
void TJ::UnpackY (complex<double>* Y, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ)
{
  int index = 0;
  for (int i = 0; i < K; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  PPsi(j, i) = Y[index]; index++;
	}
      for (int j = 0; j < J; j++)
	{
	  ZZ(j, i) = Y[index]; index++;
	}
    }
 }

// ###################################################
// Function to perform torque test on solution vectors
// ###################################################
double TJ::TorqueTest (double r, Array<complex<double>,2> YY)
{
  double torque = 0.;
  for (int i = 0; i < K; i++)
    {
      double Torque = GetTorque (r, YY, i);

      if (torque < fabs (Torque))
	torque = fabs (Torque);
    }
  
  return torque;
}

// ###############################################################################
// Function to calculate angular momentum flux associated with ith solution vector
// ###############################################################################
double TJ::GetTorque (double r, Array<complex<double>,2> YY, int i)
{
  Array<complex<double>,2> PPsi (J, K);
  Array<complex<double>,2> ZZ   (J, K);

  UnpackYY (YY, PPsi, ZZ);
  
  complex<double> I   = complex<double> (0., 1.);
  complex<double> Sum = complex<double> (0., 0.);
  double          q   = Getq (r);

  for (int j = 0; j < J; j++)
    {
      double mj  = mpol[j];
      double mnq = mj - ntor * q;
      
      Sum += (conj(ZZ(j, i)) * PPsi(j, i) - conj(PPsi(j, i)) * ZZ(j, i)) /mnq;
    }
  Sum *= I * M_PI*M_PI * ntor;

  double torque = real (Sum);

  return torque;
}

// ##################################################
// Function to calculate norms of ith solution vector
// ##################################################
void TJ::GetNorms (Array<complex<double>,2> YY, int i, double& Pnorm, double &Znorm)
{
  Array<complex<double>,2> PPsi (J, K);
  Array<complex<double>,2> ZZ   (J, K);

  UnpackYY (YY, PPsi, ZZ);
  
  double sump = 0.;
  double sumz = 0.;
  for (int j = 0; j < J; j++)
    {
      sump += real (PPsi(j, i) * conj (PPsi(j, i)));
      sumz += real (ZZ  (j, i) * conj (ZZ  (j, i)));
    }

  Pnorm = sqrt (sump);
  Znorm = sqrt (sumz);
}

