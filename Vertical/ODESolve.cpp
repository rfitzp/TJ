// ODESolve.cpp

#include "Vertical.h"

// ###################################
// Function to solve outer region odes
// ###################################
void Vertical::ODESolve ()
{
  printf ("ODE solution data:\n");
  
  // +++++++++++++++
  // Allocate memory
  // +++++++++++++++
  Array<complex<double>,2> YY (2*J, J);
  YYY.resize (2*J, J, NDIAG);
  for (int i = 0; i < 2*J; i++)
    for (int j = 0; j < J; j++)
      for (int k = 0; k < NDIAG; k++)
	YYY(i, j, k) = complex<double> (0., 0.);

  Rgrid = new double[NDIAG];
  for (int i = 0; i < NDIAG; i++)
    Rgrid[i] = EPS + (1. - EPS) * double (i) /double (NDIAG - 1);

  Pgrid = new double[NDIAG];
  for (int i = 0; i < NDIAG; i++)
    Pgrid[i] = gsl_spline_eval (Pspline, Rgrid[i], Pacc);

  hode = new double[NDIAG];
  eode = new double[NDIAG];

  Etest.resize (J, NDIAG);
  Pnorm.resize (J, NDIAG);
  Znorm.resize (J, NDIAG);

  for (int i = 0; i < NDIAG; i++)
    {
      hode[i] = h0;
      eode[i] = acc;
      for (int j = 0; j < J; j++)
	{
	  Etest (j, i) = 0.;
	  Pnorm (j, i) = 0.;
	  Znorm (j, i) = 0.;
	}
    }

  // ++++++++++++++++++++++++++++++++++++++++++++
  // Initialize solution vectors at magnetic axis
  // ++++++++++++++++++++++++++++++++++++++++++++
  double rr = EPS;
  LaunchAxis (rr, YY);

  printf ("Initializing solution vectors at magnetic axis:             r = %11.4e: Energy test = %11.4e\n", rr, EnergyTest (rr, YY));

  // +++++++++++++++++++++++++++++++++++++++++++++
  // Integrate solution vectors to plasma boundary
  // +++++++++++++++++++++++++++++++++++++++++++++
  double rx = 1.;
  int    nf = (NFIX == 0) ? 0 : NFIX;
  SegmentFixup (rr, rx, nf, YY);
  
  printf ("Integrating solution vectors to plasma boundary:            r = %11.4e: Energy test = %11.4e\n", rx, EnergyTest (rx, YY));

  int index = NDIAG - 1;

  for (int j = 0; j < 2*J; j++)
    for (int jp = 0; jp < J; jp++)
      YYY(j, jp, index) = YY(j, jp);

  hode[index] = hode[index-1];
  eode[index] = eode[index-1];
  for (int i = 0; i < J; i++)
    {
      double Pnm, Znm;
      GetNorms (YY, i, Pnm, Znm);
      
      Pnorm(i, index) = Pnm;
      Znorm(i, index) = Znm;
    }

  // ........................................................
  // Calculate energy fluxes associated with solution vectors
  // ........................................................
  for (int i = 0; i < NDIAG; i++)
    for (int jp = 0; jp < J; jp++)
      {
	complex<double> I   = complex<double> (0., 1.);
	complex<double> Sum = complex<double> (0., 0.);
	
	for (int j = 0; j < J; j++)
	  {
	    Sum += (conj(YYY(J+j, jp, i)) * YYY(j, jp, i) - conj(YYY(j, jp, i)) * YYY(J+j, jp, i));
	  }
	Sum *= I * M_PI*M_PI;

	Etest(jp, i) = real(Sum);
      }
}

// ######################################################
// Function to launch solution vectors from magnetic axis
// ######################################################
void Vertical::LaunchAxis (double r, Array<complex<double>,2> YY)
{
  Array<complex<double>,2> PPsi (J, J);
  Array<complex<double>,2> ZZ   (J, J);

  // .........................................................
  // Launch J independent solutions vectors from magnetic axis
  // .........................................................
  for (int j = 0; j < J; j++)
    {
      int    MJ  = MPOL[j];
      double mj  = mpol[j];
 
      for (int k = 0; k < J; k++)
	{
	  PPsi(k, j) = complex<double> (0., 0.);
	  ZZ  (k, j) = complex<double> (0., 0.);
	}

       if (MJ >= 0)
	 {
	   PPsi(j, j) = complex<double> (1., 0.);
	   ZZ  (j, j) = + mj;
	 }
       else if (MJ < 0)
	 {
	   PPsi(j, j) = complex<double> (1., 0.);
	   ZZ  (j, j) = - mj;
	 }
    }

  PackYY (PPsi, ZZ, YY);

  int index = 0;
 
  for (int j = 0; j < 2*J; j++)
    for (int jp = 0; jp < J; jp++)
      YYY(j, jp, index) = YY(j, jp);
  
  hode[index] = h0;
  eode[index] = acc;
 
  for (int i = 0; i < J; i++)
    {
      double Pnm, Znm;
      GetNorms (YY, i, Pnm, Znm);

      Pnorm(i, index) = Pnm;
      Znorm(i, index) = Znm;
    }
}

// ##########################################################################################################
// Function to integrate multiple solution vectors from given value of r to r = rx while performing nf fixups
// ##########################################################################################################
void Vertical::SegmentFixup (double& r, double rx, int nf, Array<complex<double>,2> YY)
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

      printf (": Energy test = %11.4e\n", EnergyTest (r, YY));
    }
}

// ##################################################################################
// Function to integrate multiple solution vectors from given value of r to to r = rx
// ##################################################################################
void Vertical::Segment (double& r, double rx, Array<complex<double>,2> YY)
{
  int              neqns = (2*J)*(J);
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
	    for (int jp = 0; jp < J; jp++)
	      YYY(j, jp, index) = YY(j, jp);
	  
	  hode[index] = h;
	  eode[index] = t_err;
	  for (int i = 0; i < J; i++)
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
void Vertical::Fixup (double r, Array<complex<double>,2> YY)
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
	}
    }
}

// ##########################################################
// Function to evaluate right-hand sides of outer region odes
// ##########################################################
void Vertical::CashKarp45Rhs (double r, complex<double>* Y, complex<double>* dYdr)
{
  Array<complex<double>,2> PPsi    (J, J);
  Array<complex<double>,2> ZZ      (J, J);
  Array<complex<double>,2> dPPsidr (J, J);
  Array<complex<double>,2> dZZdr   (J, J);
  Array<complex<double>,2> AAmmp   (J, J);
  Array<complex<double>,2> BBmmp   (J, J);
  Array<complex<double>,2> CCmmp   (J, J);
  Array<complex<double>,2> DDmmp   (J, J);

  GetMatrices (r, AAmmp, BBmmp, CCmmp, DDmmp);
  
  UnpackY (Y, PPsi, ZZ);

  for (int i = 0; i < J; i++)
    {
      for (int j = 0; j < J; j++)
	{
	  dPPsidr(j, i) = complex<double> (0., 0.);
	  dZZdr  (j, i) = complex<double> (0., 0.);

	  for (int jp = 0; jp < J; jp++)
	    {
	      dPPsidr(j, i) += (AAmmp(j, jp) * ZZ(jp, i) + BBmmp(j, jp) * PPsi(jp, i)) /r;
	      dZZdr  (j, i) += (CCmmp(j, jp) * ZZ(jp, i) + DDmmp(j, jp) * PPsi(jp, i)) /r;
	    }
	}
    }
  
  PackY (dPPsidr, dZZdr, dYdr);
}

// ###################################
// Function to pack YY solution vector
// ###################################
 void Vertical::PackYY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, Array<complex<double>,2> YY)
{
  for (int i = 0; i < J; i++)
    for (int j = 0; j < J; j++)
      {
	YY(j,   i) = PPsi(j, i);
	YY(J+j, i) = ZZ  (j, i); 
      }
}

// ######################################
// Function to unpack YY solution vectors 
// ######################################
void Vertical::UnpackYY (Array<complex<double>,2> YY, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ)
{
  for (int i = 0; i < J; i++)
    for (int j = 0; j < J; j++)
      {
	PPsi(j, i) = YY(j,   i);
	ZZ  (j, i) = YY(J+j, i);
      }
}

// ####################################
// Function to pack YY solution vectors 
// ####################################
void Vertical::PackYY (complex<double>* Y, Array<complex<double>,2> YY)
{
  int index = 0;
  for (int i = 0; i < J; i++)
    for (int j = 0; j < 2*J; j++)
      {
	YY(j, i) = Y[index]; index++;
      }
}

// #####################################
// Function to unpack YY solution vector
// #####################################
void Vertical::UnpackYY (Array<complex<double>,2> YY, complex<double>* Y)
{
  int index = 0;
  for (int i = 0; i < J; i++)
    for (int j = 0; j < 2*J; j++)
      {
	Y[index] = YY(j, i); index++;
      }
}

// ##################################
// Function to pack Y solution vector
// ##################################
void Vertical::PackY (Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ, complex<double>* Y)
{
  int index = 0;
  for (int i = 0; i < J; i++)
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
void Vertical::UnpackY (complex<double>* Y, Array<complex<double>,2> PPsi, Array<complex<double>,2> ZZ)
{
  int index = 0;
  for (int i = 0; i < J; i++)
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
// Function to perform energy test on solution vectors
// ###################################################
double Vertical::EnergyTest (double r, Array<complex<double>,2> YY)
{
  double energy = 0.;
  for (int i = 0; i < J; i++)
    {
      double Energy = GetEnergy (r, YY, i);

      if (energy < fabs (Energy))
	energy = fabs (Energy);
    }
  
  return energy;
}

// #####################################################################
// Function to calculate energy flux associated with ith solution vector
// #####################################################################
double Vertical::GetEnergy (double r, Array<complex<double>,2> YY, int i)
{
  Array<complex<double>,2> PPsi (J, J);
  Array<complex<double>,2> ZZ   (J, J);

  UnpackYY (YY, PPsi, ZZ);
  
  complex<double> I   = complex<double> (0., 1.);
  complex<double> Sum = complex<double> (0., 0.);

  for (int j = 0; j < J; j++)
    {
      Sum += (conj(ZZ(j, i)) * PPsi(j, i) - conj(PPsi(j, i)) * ZZ(j, i));
    }
  Sum *= I * M_PI*M_PI;

  double energy = real (Sum);

  return energy;
}

// ##################################################
// Function to calculate norms of ith solution vector
// ##################################################
void Vertical::GetNorms (Array<complex<double>,2> YY, int i, double& Pnorm, double &Znorm)
{
  Array<complex<double>,2> PPsi (J, J);
  Array<complex<double>,2> ZZ   (J, J);

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

