// Toroidal.cpp

#include "Vertical.h"

// ################################################################
// Function to return associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 124, Eq. (16)
// ################################################################
double Vertical::ToroidalP (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("Vertical::ToroidalP: Error - z < 1.\n");
      exit (1);
    }	

  // Calculate factor multiplying hypergeometric series
  double x  = (z - 1.) / (z + 1.);
  double d  = pow (0.5 + z/2., double (abs (n)) - 0.5);
  d        *= pow (x, double (abs (m)) /2.);
  
  if (m > 0)
    for (int j = 1; j <= m; j++)
      d *= (double (n*n - j*j + j) - 0.25) / double (j);
  else if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (j);
    }

  // Sum hypergeometric series
  double a, c, b, e, f, q, r;
  a = 0.5 - double (abs (n));
  c = 1.  + double (abs (m));
  b = a + c - 1.;
  e = 1.;
  f = 1.;
  q = 0.;
  
  do 
    {
      e *= x * (a + q) * (b + q) / (c + q) / (1. + q);
      f += e;
      q += 1.;
      
      r  = (a + q) / (a + q - 1.);
      r *= (b + q) / (b + q - 1.);
      r *= (c + q - 1.) / (c + q);
      r *= q / (q + 1.);
    }
  while (fabs (e) > 1.e-15 || fabs (r) > 1./x || q < - double (n));
  
  return f * d;
} 
  
// ################################################################
// Function to return associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// Routine sums hypergeometric series representation of function to 
// accuracy of about 1 in 10^12.
//
// Reference: Bateman Manuscript Project, Vol. I, p. 134, Eq. (41)
// ################################################################
double Vertical::ToroidalQ (int m, int n, double z)
{
  // Check argument
  if (z < 1.)
    {
      printf ("Vertical::ToroidalQ: Error - z < 1.\n");
      exit (1);
    }

  // Calculate factor multiplying hypergeometric series
  double y  = sqrt ((z*z - 1.) /z/z);
  double d  = M_PI / sqrt (2.*z); 
  int    na = abs (n);
  if (m < 0)
    {
      int mp = - m;
      for (int j = 1; j <= mp; j++)
	d /= double (n*n - j*j + j) - 0.25;
    }
  int ma = abs (m);
  if (na > 0)
    for (int k = 1; k <= na; k++)
      d *= (double (ma + k) - 0.5) /2./z / double (k);
  if (ma > 0)
    for (int l = 1; l <= ma; l++)
      d *= - y * (double (l) - 0.5);

  // Sum hypergeometric series
  double a, b, c, e, f, q, r;
  a = 0.5 * (1.5 + double (ma + na));
  b = a - 0.5;
  c = double (na + 1);
  e = 1.;
  f = 1.;
  q = 0.;
  
  do
    {
      e *= (a + q) * (b + q) / (c + q) / (1. + q) / z / z;
      f += e;
      q += 1.;
      
      r  = (a + q) / (a + q - 1.);
      r *= (b + q) / (b + q - 1.);
      r *= (c + q - 1.) / (c + q);
      r *= q / (q + 1.);
    }
  while (e > 1.e-15 || r > z*z);
  
  return f * d;
} 
  
// #############################################################
// Function to return derivative of associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// #############################################################
double Vertical::ToroidaldPdz (int m, int n, double z)
{
  return ((double (n) - 0.5) * z * ToroidalP (m, n, z) - (double (n) - 0.5 + double (m)) * ToroidalP (m, n-1, z)) /(z*z - 1.);
}

// #############################################################
// Function to return derivative of associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// #############################################################
double Vertical::ToroidaldQdz (int m, int n, double z)
{
  return ((double (n) - 0.5) * z * ToroidalQ (m, n, z) - (double (n) - 0.5 + double (m)) * ToroidalQ (m, n-1, z)) /(z*z - 1.);
}

// #############################################################
// Function to return normalized associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// #############################################################
double Vertical::NormToroidalP (int m, int n, double z)
{
  double mm = double (m);
  int    N  = abs (n);
  double nn = fabs (double (n));
  
  double Pfac = cos (nn*M_PI) * sqrt (M_PI) * gsl_sf_gamma (nn + 0.5 - mm) * pow (epsa, nn)
                /pow (2., nn - 0.5) /gsl_sf_fact (N);

  return Pfac * ToroidalP (m, n, z);
}

// ########################################################################
// Function to return normalized derivative of associated Legendre function
//  
//   P^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// ########################################################################
double Vertical::NormToroidaldPdz (int m, int n, double z)
{
  double mm = double (m);
  int    N  = abs (n);
  double nn = fabs (double (n));
  
  double Pfac = cos (nn*M_PI) * sqrt (M_PI) * gsl_sf_gamma (nn + 0.5 - mm) * pow (epsa, nn)
                /pow (2., nn - 0.5) /gsl_sf_fact (N);

  return Pfac * ToroidaldPdz (m, n, z);
}

// #############################################################
// Function to return normalized associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// #############################################################
double Vertical::NormToroidalQ (int m, int n, double z)
{
  double mm = double (m);
  int    N  = abs (n);
  double nn = fabs (double (n));
  
  double Qfac = cos (mm*M_PI) * cos (nn*M_PI) * pow (2., nn - 0.5) * gsl_sf_fact (N)
                /sqrt (M_PI) /gsl_sf_gamma (nn + 0.5 + mm) /pow (epsa, nn);

  return Qfac * ToroidalQ (m, n, z);
}

// ########################################################################
// Function to return normalized derivative of associated Legendre function
//  
//   Q^m_(n-1/2) (z)
//
//  m ... integer
//  n ... integer
//  z ... double greater than 1.0
//
// ########################################################################
double Vertical::NormToroidaldQdz (int m, int n, double z)
{
  double mm = double (m);
  int    N  = abs (n);
  double nn = fabs (double (n));
  
  double Qfac = cos (mm*M_PI) * cos (nn*M_PI) * pow (2., nn - 0.5) * gsl_sf_fact (N)
                /sqrt (M_PI) /gsl_sf_gamma (nn + 0.5 + mm) /pow (epsa, nn);

  return Qfac * ToroidaldQdz (m, n, z);
}

// ##############################################################
// Function to return hyperbolic cosine of toroidal coordinate mu
// ##############################################################
double Vertical::GetCoshMu (double R, double Z)
{
  double d1 = sqrt ((R + 1.) * (R + 1.) + Z*Z);
  double d2 = sqrt ((R - 1.) * (R - 1.) + Z*Z);

  return 0.5 * (d1 /d2 + d2 /d1);
}

// ##########################################
// Function to return toroidal coordinate eta
// ##########################################
double Vertical::GetEta (double R, double Z)
{
  double r2 = R*R + Z*Z;
  
  if (Z >= 0.)
    return + acos ((r2 - 1.) /sqrt ((r2 - 1.) * (r2 - 1.) + 4.*Z*Z));
  else
    return - acos ((r2 - 1.) /sqrt ((r2 - 1.) * (r2 - 1.) + 4.*Z*Z));
}

