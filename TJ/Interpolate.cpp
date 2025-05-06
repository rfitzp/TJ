// Interpolate.cpp

#include "TJ.h"

// #################################################################
// Functions to return interpolated values of equilibrium quantities
// #################################################################

double TJ::GetPsiN (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Pspline, 1., Pacc);
  else
    return gsl_spline_eval (Pspline, r, Pacc);
}

double TJ::Getf (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (fspline, 1., facc);
  else
    return gsl_spline_eval (fspline, r, facc);
}

double TJ::Getpp (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (ppspline, 1., ppacc);
  else
    return gsl_spline_eval (ppspline, r, ppacc);
}

double TJ::Getppp (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (pppspline, 1., pppacc);
  else
    return gsl_spline_eval (pppspline, r, pppacc);
}

double TJ::Getq (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (qspline, 1., qacc);
  else
    return gsl_spline_eval (qspline, r, qacc);
}

double TJ::Getg2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (g2spline, 1., g2acc);
  else
    return gsl_spline_eval (g2spline, r, g2acc);
}

double TJ::Getg (double r)
{
 return 1. + epsa*epsa * Getg2 (r);
}

double TJ::Gets (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (sspline, 1., sacc);
  else
    return gsl_spline_eval (sspline, r, sacc);
}

double TJ::Gets2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (s2spline, 1., s2acc);
  else
    return gsl_spline_eval (s2spline, r, s2acc);
}

double TJ::Gets0 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (s0spline, 1., s0acc);
  else
    return gsl_spline_eval (s0spline, r, s0acc);
}

double TJ::GetS1 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S1spline, 1., S1acc);
  else
    return gsl_spline_eval (S1spline, r, S1acc);
}

double TJ::GetS2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S2spline, 1., S2acc);
  else
    return gsl_spline_eval (S2spline, r, S2acc);
}

double TJ::GetS3 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S3spline, 1., S3acc);
  else
    return gsl_spline_eval (S3spline, r, S3acc);
}

double TJ::GetS4 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S4spline, 1., S4acc);
  else
    return gsl_spline_eval (S4spline, r, S4acc);
}

double TJ::GetP1 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P1spline, 1., P1acc);
  else
    return gsl_spline_eval (P1spline, r, P1acc);
}

double TJ::GetP2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P2spline, 1., P2acc);
  else
    return gsl_spline_eval (P2spline, r, P2acc);
}

double TJ::GetP3 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P3spline, 1., P3acc);
  else
    return gsl_spline_eval (P3spline, r, P3acc);
}

double TJ::Getne (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (nespline, 1., neacc);
  else
    return gsl_spline_eval (nespline, r, neacc);
}

double TJ::GetTe (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Tespline, 1., Tepacc);
  else
    return gsl_spline_eval (Tespline, r, Tepacc);
}
 
double TJ::Getnep (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (nepspline, 1., nepacc);
  else
    return gsl_spline_eval (nepspline, r, nepacc);
}

double TJ::GetTep (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Tepspline, 1., Tepacc);
  else
    return gsl_spline_eval (Tepspline, r, Tepacc);
}
 
double TJ::GetHn (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (HHspline[n], 1., HHacc[n]);
  else
    return gsl_spline_eval (HHspline[n], r, HHacc[n]);
}

double TJ::GetHnp (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (HPspline[n], 1., HPacc[n]);
  else
    return gsl_spline_eval (HPspline[n], r, HPacc[n]);
}

double TJ::GetVn (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (VVspline[n], 1., VVacc[n]);
  else
    return gsl_spline_eval (VVspline[n], r, VVacc[n]);
}

double TJ::GetVnp (int n, double r)
{
   if (r >= 1.)
     return gsl_spline_eval (VPspline[n], 1., VPacc[n]);
   else
     return gsl_spline_eval (VPspline[n], r, VPacc[n]);
}

double TJ::GetDI (double r)
{
  double pp = Getpp (r);
  double q  = Getq (r);
  double s  = Gets (r);

  return - 0.25 - epsa*epsa * 2. * r*pp * (1. - q*q) /s/s; 
}

double TJ::GetDR (double r)
{
  double pp  = Getpp (r);
  double H1p = GetHnp (1, r);
  double q   = Getq (r);
  double s   = Gets (r);

  return - epsa*epsa * 2. * r*pp * (1. - q*q) /s/s - epsa*epsa * 2. * pp * q*q * H1p /s; 
}

double TJ::GetFlarge (double r, int m)
{
  double f   = Getf (r);
  double F   = Psi[Nr];
  double s   = Gets (r);
  double S1  = GetS1 (r);
  double H1  = GetHn (1, r);
  double DI  = GetDI (r);
  double nuL = 0.5 - sqrt (- DI);
  double nuS = 0.5 + sqrt (- DI);

  // Assume that DCON psi is PsiN
  double rho   = f /F;
  double dPNdP = 1.;

  double mm   = double (m);
  double eps2 = epsa*epsa;
  double r2   = r*r;
  double mm2  = mm*mm;
  double nt2  = ntor*ntor;
  double Lmm  = mm2 + eps2 * mm2 * (- 0.75*r2 + H1 + S1) + eps2*nt2*r2;

  return pow (rho, nuL - 1.) * sqrt ((nuS - nuL) /Lmm) * s * mm * F * dPNdP;
}

double TJ::GetFsmall (double r, int m)
{
  double f   = Getf (r);
  double F   = Psi[Nr];
  double s   = Gets (r);
  double S1  = GetS1 (r);
  double H1  = GetHn (1, r);
  double DI  = GetDI (r);
  double nuL = 0.5 - sqrt (- DI);
  double nuS = 0.5 + sqrt (- DI);

  // Assume that DCON psi is PsiN
  double rho   = f /F;
  double dPNdP = 1.;

  double mm   = double (m);
  double eps2 = epsa*epsa;
  double r2   = r*r;
  double mm2  = mm*mm;
  double nt2  = ntor*ntor;
  double Lmm  = mm2 + eps2 * mm2 * (- 0.75*r2 + H1 + S1) + eps2*nt2*r2;

  return pow (rho, nuS - 1.) * sqrt ((nuS - nuL) /Lmm) * s * mm * F * dPNdP;
}
