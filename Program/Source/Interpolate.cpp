// Interpolate.cpp

#include "TJ.h"

// #################################################################
// Functions to return interpolated values of equilibrium quantities
// #################################################################

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

double TJ::GetS1 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S1spline, 1., S1acc);
  else
    return gsl_spline_eval (S1spline, r, S1acc);
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
