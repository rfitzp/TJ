// Resonant.cpp

#include "TJ.h"

// #################################################################
// Function to calculate coefficients L1, P1, and T1 at given radius
// #################################################################
void TJ::GetL1P1T1 (double rm, double mm, int km, double& sm, double& L1, double& P1, double& T1)
{
  double qm, s2m, qp, qpp;

  qm  = Getq  (rm);
  sm  = Gets  (rm);
  s2m = Gets2 (rm);
  qp  = qm * sm  /rm;
  qpp = qm * s2m /rm/rm;
  T1  = 1. + 0.5 * rm * qpp /qp;

  complex<double> Lmmp, Mmmp, Nmmp, Pmmp;
  
  GetMatrices (rm, km, km, Lmmp, Mmmp, Nmmp, Pmmp);
  double LM = real (Lmmp);
  double PM = real (Pmmp);

  GetMatrices (rm + EPSF, km, km, Lmmp, Mmmp, Nmmp, Pmmp);
  double LMP = real (Lmmp);
  double PMP = real (Pmmp);

  GetMatrices (rm - EPSF, km, km, Lmmp, Mmmp, Nmmp, Pmmp);
  LMP -= real (Lmmp);
  PMP -= real (Pmmp);
  
  LMP /= 2.*EPSF;
  PMP /= 2.*EPSF;
  
  L1 = (rm /mm /sm) * (0.5 * (qpp /qp) * LM - LMP);
  P1 = (rm /mm /sm) * (0.5 * (qpp /qp) * PM - PMP);
}

// #########################################################################
// Function to calculate coefficients L1k, M1k, N1k, and P1k at given radius
// #########################################################################
void TJ::GetL1kM1kN1kP1k (double rm, double mm, int km, complex<double>* L1k, complex<double>* M1k,
			  complex<double>* N1k, complex<double>* P1k)
  {
    double qm, sm, s2m, qp, qpp;

    qm  = Getq  (rm);
    sm  = Gets  (rm);
    s2m = Gets2 (rm);
    qp  = qm * sm  /rm;
    qpp = qm * s2m /rm/rm;

    complex<double> Lmmp, Mmmp, Nmmp, Pmmp, LM, MM, NM, PM, LMP, MMP, NMP, PMP; 
    for (int k = 0; k < J; k++)
      {
	GetMatrices (rm, k, km, Lmmp, Mmmp, Nmmp, Pmmp);
	LM = Lmmp;
	MM = Mmmp;
	NM = Nmmp;
	PM = Pmmp;

	GetMatrices (rm + EPSF, k, km, Lmmp, Mmmp, Nmmp, Pmmp);
	LMP = Lmmp;
	MMP = Mmmp;
	NMP = Nmmp;
	PMP = Pmmp;

	GetMatrices (rm - EPSF, k, km, Lmmp, Mmmp, Nmmp, Pmmp);
	LMP -= Lmmp;
	MMP -= Mmmp;
	NMP -= Nmmp;
	PMP -= Pmmp;
	
	LMP /= 2.*EPSF;
	MMP /= 2.*EPSF;
	NMP /= 2.*EPSF;
	PMP /= 2.*EPSF;

	L1k[k] = (rm /mm /sm) * (0.5 * (qpp/qp) * LM - LMP);
	M1k[k] = (rm /mm /sm) * (0.5 * (qpp/qp) * MM - MMP);
	N1k[k] = (rm /mm /sm) * (0.5 * (qpp/qp) * NM - NMP);
	P1k[k] = (rm /mm /sm) * (0.5 * (qpp/qp) * PM - PMP);
      }
}

// ##############################################################
// Function to jump solution vectors across rational surface jres
// ##############################################################
void TJ::JumpRational (int jres, double& r, Array<complex<double>,2> YY)
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

  // ............................
  // Calculate L1k, M1k, N1k, P1k
  // ............................
  complex<double>* L1k = new complex<double>[J];
  complex<double>* M1k = new complex<double>[J];
  complex<double>* N1k = new complex<double>[J];
  complex<double>* P1k = new complex<double>[J];
  GetL1kM1kN1kP1k (rm, mm, km, L1k, M1k, N1k, P1k);

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

  // ++++++++++++++++++++
  // Finite pressure case
  // ++++++++++++++++++++
  if (fabs (nuL) > NULC)
    {
      double dnuL = pow (DEL, nuL);
      double dnuS = pow (DEL, nuS);

      // ................
      // Calculate bL, bS
      // ................
      double bL = nuL /L0;
      double bS = nuS /L0;
      
      // ............
      // Perform sums
      // ............
      complex<double> Sum1 = complex<double> (0., 0.);
      complex<double> Sum2 = complex<double> (0., 0.);

      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      double kk = mpol[k] - mpol[km];
	      
	      Sum1 += (Lmat(km, k) * Pmat(k, km) + Pmat(km, k) * Lmat(k, km) + Mmat(km, k) * Mmat(k, km) + Nmat(km, k) * Nmat(k, km)
		       + bL * (Lmat(km, k) * Nmat(k, km) + Mmat(km, k) * Lmat(k, km))
		       + (Nmat(km, k) * Pmat(k, km) + Pmat(km, k) * Mmat(k, km)) /bL) /kk;
	      
	      Sum2 += (  (nuL + 1.) * (Pmat(km, k) * Lmat(k, km) + Nmat(km, k) * Nmat(k, km))
		       + (nuL - 1.) * (Lmat(km, k) * Pmat(k, km) + Mmat(km, k) * Mmat(k, km))
		       + bL * (nuL - 1.) * (Lmat(km, k) * Nmat(k, km) + Mmat(km, k) * Lmat(k, km))
		       + (nuL + 1.) * (Nmat(km, k) * Pmat(k, km) + Pmat(km, k) * Mmat(k, km)) /bL) /kk;
	    }
	}
      Sum1 /= - 2. * mm * sm * rm * nuL;
      Sum2 /= - 2. * mm * sm * rm * nuL * L0;
      
      // ....................
      // Calculate lamL, gamL
      // ....................
      complex<double> lamL = (L0 * P1 /nuL + T1 + nuL * (L1 /L0 - 2.)) /2./rm + Sum1;
      complex<double> gamL = ((1. + nuL) * (P1 /nuL + T1 /L0 - nuL /L0) + P0 * (L1 /L0 - 1.)) /2./rm + Sum2;

      // ................
      // Calculate ak, bk
      // ................
      complex<double>* ak = new complex<double>[J];
      complex<double>* bk = new complex<double>[J];
      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      ak[k] = - (Lmat(k, km) /L0  + Mmat(k, km) /nuL) /mm /sm;
	      bk[k] = - (Pmat(k, km) /nuL + Nmat(k, km) /L0 ) /mm /sm;
	    }
	  else
	    {
	      ak[k] = complex<double> (0., 0.);
	      bk[k] = complex<double> (0., 0.);
	    }
	}

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

      // ................
      // Calculate ck, dk
      // ................
      complex<double>* ck = new complex<double>[J];
      complex<double>* dk = new complex<double>[J];

      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      ck[k] = - nuL * ak[k] + L1k[k] * bL + M1k[k]
		- (rm /mm /sm) * (Lmat(k, km) * gamL + Mmat(k, km) * lamL);
	      dk[k] = - (mm*sm /(mpol[k] - mm) + nuL) * bk[k] + N1k[k] * bL + P1k[k]
		- (rm /mm /sm) * (Nmat(k, km) * gamL + Pmat(k, km) * lamL);
	      
	      for (int kk = 0; kk < J; kk++)
		{
		  if (kk != km)
		    {
		      double kx = mpol[kk] - mm;
		      
		      ck[k] += (Lmat(k, kk) * bk[kk] + Mmat(k, kk) * ak[kk]) /kx;
		      dk[k] += (Nmat(k, kk) * bk[kk] + Pmat(k, kk) * ak[kk]) /kx;
		    }
		}
	      
	      ck[k] /= rm * (1. + nuL);
	      dk[k] /= rm * (1. + nuL);
	    }
	  else
	    {
	      ck[k] = complex<double> (0., 0.);
	      dk[k] = complex<double> (0., 0.);
	    }
	}
 
      // ..................................................
      // Calculate AL, AS, AC, BC, Pbar, Zbar, Pbar1, Zbar1
      // ..................................................
      double                   Tminus, Tplus;
      complex<double>          AL, AS, AC, BC; 
      Array<complex<double>,2> PPsi (J, K), ZZ (J, K);
      complex<double>*         Pbar  = new complex<double>[J];
      complex<double>*         Zbar  = new complex<double>[J];
      complex<double>*         Pbar1 = new complex<double>[J];
      complex<double>*         Zbar1 = new complex<double>[J];
	  
      UnpackYY (YY, PPsi, ZZ);

      for (int j = 0; j < K; j++)
	{
	  AL = PPsi(km, j) /dnuL;
	  AS = complex<double> (0., 0.);

	  for (int k = 0; k < J; k++)
	    {
	      Pbar1[k] = complex<double> (0., 0.);
	      Zbar1[k] = complex<double> (0., 0.);
	    }

	  Tminus = GetTorque (rm - DEL, PPsi, ZZ, j);

	  for (int iter = 0; iter < ITERMAX; iter++)
	    {
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      Pbar[k] = PPsi(k, j) - (ak[k] - ck[k]*DEL) * AL * dnuL + akt[k] * AS * dnuS + Pbar1[k] * DEL;
		      Zbar[k] = ZZ  (k, j) - (bk[k] - dk[k]*DEL) * AL * dnuL + bkt[k] * AS * dnuS + Zbar1[k] * DEL;
		    }
		  else
		    {
		      Pbar[k] = complex<double> (0., 0.);
		      Zbar[k] = complex<double> (0., 0.);
		    }
		}
	      
	      AC = complex<double> (0., 0.);
	      BC = complex<double> (0., 0.);
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      double kk = mpol[k] - mpol[km];
		      
		      AC += (Nmat(km, k) * Zbar[k] + Pmat(km, k) * Pbar[k]) /kk;
		      BC += (Lmat(km, k) * Zbar[k] + Mmat(km, k) * Pbar[k]) /kk;
		    }
		}
	      AC /= - rm*P0;
	      BC /= - rm*L0;
	      BC +=   AC/L0;

	      // ......................
	      // Calculate Pbar1, Zbar1
	      // ......................
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      Pbar1[k] = - (rm /mm /sm) * (Lmat(k, km) * BC + Mmat(k, km) * AC);
		      Zbar1[k] = - (mm * sm /(mpol[k] - mm)) * Zbar[k]
			- (rm /mm /sm) * (Nmat(k, km) * BC + Pmat(k, km) * AC);
		      
		      for (int kk = 0; kk < J; kk++)
			{
			  if (kk != km)
			    {
			      double kx = mpol[kk] - mm;
			      
			      Pbar1[k] += (Lmat(k, kk) * Zbar[kk] + Mmat(k, kk) * Pbar[kk]) /kx;
			      Zbar1[k] += (Nmat(k, kk) * Zbar[kk] + Pmat(k, kk) * Pbar[kk]) /kx;
			    }
			}
		      
		      Pbar1[k] /= rm;
		      Zbar1[k] /= rm;
		    }
		  else
		    {
		      Pbar1[k] = complex<double> (0., 0.);
		      Zbar1[k] = complex<double> (0., 0.);
		    }
		}  

	      AS = - (ZZ(km, j) - bL * PPsi(km, j) + DEL*(BC - bL*AC) + DEL*(gamL - bL*lamL)*AL*dnuL)
		/(bS - bL) /dnuS;
	      AL =   (PPsi(km, j) + AS*dnuS + AC*DEL) /(1. - DEL*lamL) /dnuL;
	    }

	  Pi(jres, j) = pow (rm, nuL) * sqrt ((nuS - nuL) /real(Lmat(km, km))) * AL;

	  // ............
	  // Perform jump
	  // ............
	  for (int k = 0; k < J; k++)
	    {
	      if (k == km)
		{
		  PPsi(k, j) += 2. * DEL * (AL*dnuL * lamL + AC) + 2.      * AS * dnuS;
		  ZZ  (k, j) += 2. * DEL * (AL*dnuL * gamL + BC) + 2. * bS * AS * dnuS;
		}
	      else
		{
		  PPsi(k, j) += 2. * DEL * (AL * dnuL * ck[k] + Pbar1[k]) + 2. * akt[k] * AS * dnuS;
		  ZZ  (k, j) += 2. * DEL * (AL * dnuL * dk[k] + Zbar1[k]) + 2. * bkt[k] * AS * dnuS;
		} 
	    }

	  Tplus = GetTorque (rm + DEL, PPsi, ZZ, j);

	  if (j < J)
	    printf ("Surface %2d nuL = %10.3e m = %3d AL = (%10.3e, %10.3e) AS = (%10.3e, %10.3e) Tminus = %10.3e Tplus = %10.3e\n",
		    jres+1, nuL, MPOL[j], real(AL), imag(AL), real(AS), imag(AS), Tminus, Tplus);
	  else
	    printf ("Surface %2d nuL = %10.3e m = %3d AL = (%10.3e, %10.3e) AS = (%10.3e, %10.3e) Tminus = %10.3e Tplus = %10.3e\n",
		    jres+1, nuL, mres[j-J], real(AL), imag(AL), real(AS), imag(AS), Tminus, Tplus);
	}

      // .............
      // Jump r and YY
      // .............
      r += 2.*DEL;
      PackYY (PPsi, ZZ, YY);

      delete[] Pbar; delete[] Zbar; delete[] ak;    delete[] bk; delete[] akt; delete[] bkt;
      delete[] L1k;  delete[] M1k;  delete[] N1k;   delete[] P1k;
      delete[] ck;   delete[] dk;   delete[] Pbar1; delete[] Zbar1;
     }
 else
    // +++++++++++++++++
    // Low pressure case
    // +++++++++++++++++
    {
      // ........................
      // Calculate bS, gamL, epsL
      // ........................
      double bL   = nuL /L0;
      double bS   = nuS /L0;
      double gamL = P1 * (1. + nuL) /rm + nuL * T1 /rm /L0;
      double xiL  = nuL * (L1 /L0 - 1.) /rm;
      double logd = log (DEL);

      // ............
      // Perform sums
      // ............
      complex<double> Sum1 = complex<double> (0., 0.);
      complex<double> Sum2 = complex<double> (0., 0.);
      complex<double> Sum3 = complex<double> (0., 0.);
      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      double kk = mpol[k] - mpol[km];
	      
	      Sum1 += (Lmat(km, k) * Pmat(k, km) + Mmat(km, k) * Mmat(k, km)) /kk;
	      Sum2 += (Lmat(km, k) * Nmat(k, km) + Mmat(km, k) * Lmat(k, km)) /kk;
	      Sum3 += (Nmat(km, k) * Pmat(k, km) + Pmat(km, k) * Mmat(k, km)) /kk;
	    }
	}
      Sum1 /= - mm * sm * rm;
      Sum2 /= - mm * sm * rm * L0;
      Sum3 /= - mm * sm * rm * 2. /L0;

      // .........................
      // Calculate lamL, muL, delL
      // .........................
      complex<double> lamL = P1 * L0 * (1. + nuL) /rm + nuL * T1 /rm + (Sum1 + nuL * Sum2);
      complex<double> muL  = Sum3;
      complex<double> delL = muL /L0;
 
      // ................
      // Calculate ak, bk
      // ................
      complex<double>* ak = new complex<double>[J];
      complex<double>* bk = new complex<double>[J];
      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      ak[k] = - (nuL * Lmat(k, km) /L0 + Mmat(k, km)) /mm /sm;
	      bk[k] = - (nuL * Nmat(k, km) /L0 + Pmat(k, km)) /mm /sm; 
	    }
	  else
	    {
	      ak[k] = complex<double> (0., 0.);
	      bk[k] = complex<double> (0., 0.);
	    }
	}

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
	      bkt[k] = - (Nmat(k, km) /L0  + Pmat(k, km) /nuS) /mm /sm;
	    }
	  else
	    {
	      akt[k] = complex<double> (0., 0.);
	      bkt[k] = complex<double> (0., 0.);
	    }
	}

      // ....................................
      // Calculate ck, c1k, c2k, dk, d1k, d2k
      // ....................................
      complex<double>* ck  = new complex<double>[J];
      complex<double>* c1k = new complex<double>[J];
      complex<double>* c2k = new complex<double>[J];
      complex<double>* dk  = new complex<double>[J];
      complex<double>* d1k = new complex<double>[J];
      complex<double>* d2k = new complex<double>[J];

      for (int k = 0; k < J; k++)
	{
	  if (k != km)
	    {
	      ck [k] = - ak[k] + L1k[k] * bL + M1k[k] * (1. - nuL)
		+ (rm /mm /sm) * (Lmat(k, km) * (gamL - 2.*delL) + Mmat(k, km) * (2.*lamL - 6.*muL - xiL));
	      c1k[k] = M1k[k] * nuL
		- (rm /mm /sm) * (Lmat(k, km) * (gamL - 2.*delL) + Mmat(k, km) * (lamL - 4.*muL));
	      c2k[k] = - (rm /mm /sm) * (Lmat(k, km) * delL + Mmat(k, km) * muL);
	      
	      dk [k] = - (1. - mm*sm /(mpol[k] - mm)) * bk[k] + N1k[k] * bL + P1k[k] * (1. - nuL)
			  + (rm /mm /sm) * (Nmat(k, km) * (gamL - 2.*delL) + Pmat(k, km) * (2.*lamL - 6.*muL - xiL));
	      d1k[k] = - (mm*sm /(mpol[k] - mm)) * bk[k] + P1k[k] * nuL
		- (rm /mm /sm) * (Nmat(k, km) * (gamL - 2.*delL) + Pmat(k, km) * (lamL - 4.*muL));
	      d2k[k] = - (rm /mm /sm) * (Nmat(k, km) * delL + Pmat(k, km) * muL);
			  
	      for (int kk = 0; kk < J; kk++)
		{
		  if (kk != km)
		    {
		      double kx = mpol[kk] - mm;
		      
		      ck [k] += - (Lmat(k, kk) * bk[kk] + Mmat(k, kk) * ak[kk]) /kx;
		      c1k[k] +=   (Lmat(k, kk) * bk[kk] + Mmat(k, kk) * ak[kk]) /kx;

		      dk [k] += - (Nmat(k, kk) * bk[kk] + Pmat(k, kk) * ak[kk]) /kx;
		      d1k[k] +=   (Nmat(k, kk) * bk[kk] + Pmat(k, kk) * ak[kk]) /kx;
		    }
		}
	      
	      ck [k] /= rm;
	      c1k[k] /= rm;
	      c2k[k] /= rm;

	      dk [k] /= rm;
	      d1k[k] /= rm;
	      d2k[k] /= rm;
	    }
	  else
	    {
	      ck [k] = complex<double> (0., 0.);
	      c1k[k] = complex<double> (0., 0.);
	      c2k[k] = complex<double> (0., 0.);

	      dk [k] = complex<double> (0., 0.);
	      d1k[k] = complex<double> (0., 0.);
	      d2k[k] = complex<double> (0., 0.);
	    }
	}
      
      // ....................................................................
      // Calculate AL, AS, AC, AD, BD, Pbar, Zbar, Pbar1, Zbar1, Pbar2, Zbar2
      // ....................................................................
      complex<double>          AL, AS, AC, AD, BD; 
      Array<complex<double>,2> PPsi(J, K), ZZ(J, K);
      complex<double>*         Pbar  = new complex<double>[J];
      complex<double>*         Zbar  = new complex<double>[J];
      complex<double>*         Pbar1 = new complex<double>[J];
      complex<double>*         Zbar1 = new complex<double>[J];
      complex<double>*         Pbar2 = new complex<double>[J];
      complex<double>*         Zbar2 = new complex<double>[J];
      double                   Tminus, Tplus;

      UnpackYY (YY, PPsi, ZZ);

      for (int j = 0; j < K; j++)
	{
	  AL = PPsi(km, j) / (1. + nuL * logd);
	  AS = complex<double> (0., 0.);

	  for (int k = 0; k < J; k++)
	    {
	      Pbar1[k] = complex<double> (0., 0.);
	      Zbar1[k] = complex<double> (0., 0.);
	      Pbar2[k] = complex<double> (0., 0.);
	      Zbar2[k] = complex<double> (0., 0.);
	    }

	  Tminus = GetTorque (rm - DEL, PPsi, ZZ, j);

	  for (int iter = 0; iter < ITERMAX; iter++)
	    {
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      Pbar[k] = PPsi(k, j) - (ak[k]*logd - DEL * (ck[k] + c1k[k]*logd + c2k[k]*logd*logd)) * AL
			+ akt[k] * AS * DEL + (Pbar1[k] + Pbar2[k]*logd) * DEL;
		      Zbar[k] = ZZ  (k, j) - (bk[k]*logd - DEL * (dk[k] + d1k[k]*logd + d2k[k]*logd*logd)) * AL
			+ bkt[k] * AS * DEL + (Zbar1[k] + Zbar2[k]*logd) * DEL;
		    }
		  else
		    {
		      Pbar[k] = complex<double> (0., 0.);
		      Zbar[k] = complex<double> (0., 0.);
		    }
		}
	      
	      AC = complex<double> (0., 0.);
	      AD = complex<double> (0., 0.);
	      BD = complex<double> (0., 0.);
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      double kk = mpol[k] - mpol[km];
		      
		      AC += (Lmat(km, k) * Zbar[k] + Mmat(km, k) * Pbar[k]) /kk;
		      AD += (Nmat(km, k) * Zbar[k] + Pmat(km, k) * Pbar[k]) /kk;
		    }
		}
	      AC /= rm;
	      AD /= rm /L0;
	      AD += - nuL * AC;
	      BD  = AD /L0;

	      // ....................................
	      // Calculate Pbar1, Zbar1, Pbar2, Zbar2
	      // ....................................
	      for (int k = 0; k < J; k++)
		{
		  if (k != km)
		    {
		      Pbar1[k] = - (rm /mm /sm) * (- Lmat(k, km) * BD + Mmat(k, km) * (AC - 2.*AD));
		      Zbar1[k] = - (mm * sm /(mpol[k] - mm)) * Zbar[k]
			- (rm /mm /sm) * (- Nmat(k, km) * BD + Pmat(k, km) * (AC - 2.*AD));

		      Pbar2[k] = - (rm /mm /sm) * (Lmat(k, km) * BD + Mmat(k, km) * AD);
		      Zbar2[k] = - (rm /mm /sm) * (Nmat(k, km) * BD + Pmat(k, km) * AD);

		       for (int kk = 0; kk < J; kk++)
			{
			  if (kk != km)
			    {
			      double kx = mpol[kk] - mm;
			      
			      Pbar1[k] += (Lmat(k, kk) * Zbar[kk] + Mmat(k, kk) * Pbar[kk]) /kx;
			      Zbar1[k] += (Nmat(k, kk) * Zbar[kk] + Pmat(k, kk) * Pbar[kk]) /kx;
			    }
			}

		       Pbar1[k] /= rm;
		       Zbar1[k] /= rm;

		       Pbar2[k] /= rm;
		       Zbar2[k] /= rm;
		    }
		  else
		    {
		      Pbar1[k] = complex<double> (0., 0.);
		      Zbar1[k] = complex<double> (0., 0.);

		      Pbar2[k] = complex<double> (0., 0.);
		      Zbar2[k] = complex<double> (0., 0.);
		    }
		}
	      
	      AS = - (ZZ(km, j) - AL * bL + DEL*logd*(BD + (gamL + delL*logd)*AL)) /bS /DEL;
	      AL =   (PPsi(km, j) + DEL*(AS + AC + AD*(logd - 1.)))
		/(1. + nuL * logd - DEL*(lamL*(logd - 1.) + muL*(logd*logd - 2.*logd + 2.) + xiL));
	    }
	  
	  Pi(jres, j) = sqrt (1. /real(Lmat(km, km))) * AL;
	  
	  // ............
	  // Perform jump
	  // ............
	  for (int k = 0; k < J; k++)
	    {
	      if (k == km)
		{
		  PPsi(k, j) += 2. * DEL * (AL * (xiL + lamL * (logd - 1.) + muL*(logd*logd - 2.*logd + 2.)) + AS + AC + AD * (logd - 1.));
		  ZZ  (k, j) += 2. * DEL * (AL * logd*(gamL + delL*logd) + AS * bS + BD * logd);
		}
	      else
		{
		  PPsi(k, j) += 2. * DEL * (AL * (ck[k] + c1k[k]*logd + c2k[k]*logd*logd) + akt[k] * AS + Pbar1[k] + Pbar2[k]*logd);
		  ZZ  (k, j) += 2. * DEL * (AL * (dk[k] + d1k[k]*logd + d2k[k]*logd*logd) + bkt[k] * AS + Zbar1[k] + Zbar2[k]*logd);
		} 
	    }

	  Tplus = GetTorque (rm + DEL, PPsi, ZZ, j);

	  if (j < J)
	    printf ("Surface %2d nuL = %10.3e m %3d AL = (%10.3e, %10.3e) AS = (%10.3e, %10.3e) Tminus = %10.3e Tplus = %10.3e\n",
		    jres+1, nuL, MPOL[j], real(AL), imag(AL), real(AS), imag(AS), Tminus, Tplus);
	  else
	     printf ("Surface %2d nuL = %10.3e m %3d AL = (%10.3e, %10.3e) AS = (%10.3e, %10.3e) Tminus = %10.3e Tplus = %10.3e\n",
		     jres+1, nuL, mres[j-J], real(AL), imag(AL), real(AS), imag(AS), Tminus, Tplus);
	}

      // .............
      // Jump r and YY
      // .............
      r += 2.*DEL;
      PackYY (PPsi, ZZ, YY);
         	  
      delete[] Pbar;  delete[] Zbar;  delete[] ak;    delete[] bk;   delete[] akt; delete[] bkt;
      delete[] ck;    delete[] c1k;   delete[] c2k;   delete[] dk;   delete[] d1k, delete[] d2k;
      delete[] Pbar1; delete[] Zbar1; delete[] Pbar2; delete[] Zbar2;
    }
}

// ###############################################################################
// Function to calculate angular momentum flux associated with ith solution vector
// ###############################################################################
double TJ::GetTorque (double r, Array<complex<double>,2> Psi, Array<complex<double>,2> Z, int i)
{
  complex<double> I   = complex<double> (0., 1.);
  complex<double> Sum = complex<double> (0., 0.);
  double          q   = Getq (r);

  for (int j = 0; j < J; j++)
    {
      double mj  = mpol[j];
      double mnq = mj - ntor * q;
      
      Sum += (conj(Z(j, i)) * Psi(j, i) - conj(Psi(j, i)) * Z(j, i)) /mnq;
    }
  Sum *= I * M_PI*M_PI * ntor;

  double torque = real (Sum);

  return torque;
}

