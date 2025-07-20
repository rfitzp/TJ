// Interpolate.cpp

#include "Vertical.h"

// #################################################################
// Functions to return interpolated values of equilibrium quantities
// #################################################################

double Vertical::GetPsiN (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Pspline, 1., Pacc);
  else
    return gsl_spline_eval (Pspline, r, Pacc);
}

double Vertical::Getf (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (fspline, 1., facc);
  else
    return gsl_spline_eval (fspline, r, facc);
}

double Vertical::Getp (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (p2spline, 1., ppacc);
  else
    return gsl_spline_eval (p2spline, r, ppacc);
}

double Vertical::Getpp (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (ppspline, 1., ppacc);
  else
    return gsl_spline_eval (ppspline, r, ppacc);
}

double Vertical::Getppp (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (pppspline, 1., pppacc);
  else
    return gsl_spline_eval (pppspline, r, pppacc);
}

double Vertical::Getq (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (qspline, 1., qacc);
  else
    return gsl_spline_eval (qspline, r, qacc);
}

double Vertical::Getg2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (g2spline, 1., g2acc);
  else
    return gsl_spline_eval (g2spline, r, g2acc);
}

double Vertical::Getg (double r)
{
 return 1. + epsa*epsa * Getg2 (r);
}

double Vertical::Gets (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (sspline, 1., sacc);
  else
    return gsl_spline_eval (sspline, r, sacc);
}

double Vertical::Gets2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (s2spline, 1., s2acc);
  else
    return gsl_spline_eval (s2spline, r, s2acc);
}

double Vertical::Gets0 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (s0spline, 1., s0acc);
  else
    return gsl_spline_eval (s0spline, r, s0acc);
}

double Vertical::GetS1 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S1spline, 1., S1acc);
  else
    return gsl_spline_eval (S1spline, r, S1acc);
}

double Vertical::GetS2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S2spline, 1., S2acc);
  else
    return gsl_spline_eval (S2spline, r, S2acc);
}

double Vertical::GetS3 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S3spline, 1., S3acc);
  else
    return gsl_spline_eval (S3spline, r, S3acc);
}

double Vertical::GetS4 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (S4spline, 1., S4acc);
  else
    return gsl_spline_eval (S4spline, r, S4acc);
}

double Vertical::GetP1 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P1spline, 1., P1acc);
  else
    return gsl_spline_eval (P1spline, r, P1acc);
}

double Vertical::GetP2 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P2spline, 1., P2acc);
  else
    return gsl_spline_eval (P2spline, r, P2acc);
}

double Vertical::GetP3 (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (P3spline, 1., P3acc);
  else
    return gsl_spline_eval (P3spline, r, P3acc);
}

double Vertical::Getne (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (nespline, 1., neacc);
  else
    return gsl_spline_eval (nespline, r, neacc);
}

double Vertical::GetTe (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Tespline, 1., Tepacc);
  else
    return gsl_spline_eval (Tespline, r, Tepacc);
}
 
double Vertical::Getnep (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (nepspline, 1., nepacc);
  else
    return gsl_spline_eval (nepspline, r, nepacc);
}

double Vertical::GetTep (double r)
{
  if (r >= 1.)
    return gsl_spline_eval (Tepspline, 1., Tepacc);
  else
    return gsl_spline_eval (Tepspline, r, Tepacc);
}
 
double Vertical::GetHn (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (HHspline[n], 1., HHacc[n]);
  else
    return gsl_spline_eval (HHspline[n], r, HHacc[n]);
}

double Vertical::GetHnp (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (HPspline[n], 1., HPacc[n]);
  else
    return gsl_spline_eval (HPspline[n], r, HPacc[n]);
}

double Vertical::GetVn (int n, double r)
{
  if (r >= 1.)
    return gsl_spline_eval (VVspline[n], 1., VVacc[n]);
  else
    return gsl_spline_eval (VVspline[n], r, VVacc[n]);
}

double Vertical::GetVnp (int n, double r)
{
   if (r >= 1.)
     return gsl_spline_eval (VPspline[n], 1., VPacc[n]);
   else
     return gsl_spline_eval (VPspline[n], r, VPacc[n]);
}


