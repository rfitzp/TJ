// Interpolate.cpp

// PROGRAM ORGANIZATION:
//
// double Flux:: Interpolate                  (int I, double* X, double* Y, double x, int order)
// double Flux:: InterpolateCubic             (double* X, double* Y, double x, int i0, int i1, int i2, int order)
// double Flux:: InterpolateQuartic           (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
// double Flux:: Interpolate1                 (int I, double* X, gsl_matrix * Y, int i, double x, int order)
// double Flux:: InterpolateCubic1            (double* X, gsl_matrix* Y, int i, double x, int i0, int i1, int i2, int order)
// double Flux:: InterpolateQuartic1          (double* X, gsl_matrix * Y, int i, double x, int i0, int i1, int i2, int i3, int order)
// double Flux:: Interpolate                  (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2>, int order, int cubic)
// double Flux:: InterpolateCubicCubic        (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2>, int i0, int i1, int i2, int j0, int j1, int j2, int order)
// double Flux:: InterpolateQuarticCubic      (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2>, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int order)
// double Flux:: InterpolateCubicQuartic      (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2>, int i0, int i1, int i2, int j0, int j1, int j2, int j3, int order)
// double Flux:: InterpolateQuarticQuartic    (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2>, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order)
// double Flux:: GetPsi                       (double r, double z)
// double Flux:: GetPsiR                      (double r, double z)
// double Flux:: GetPsiZ                      (double r, double z)

#include "Flux.h"

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// order = 2: d^2Y/dx^2
// ########################################################
double Flux::Interpolate (int I, double* X, double* Y, double x, int order)
{
  int index, cntrl = 0;
  
  if (x < X[0])
    {
      index = 0;
      cntrl = 2;
    }
  else if (x >= X[I-1])
    {
      index = I - 2;
      cntrl = 3;
    }
  else
    {
      for (int i = 0; i < I-1; i++)
	if (x >= X[i] && x < X[i+1])
	  {
	    index = i;

	    if (index == 0)
	      cntrl = 2;
	    else if (index == I-2)
	      cntrl = 3;
	    else
	      cntrl = 1;

	    break;
	  }
    }

  double val;
  if (cntrl == 1)
    val = InterpolateQuartic (X, Y, x, index-1, index,   index+1, index+2, order);
  else if (cntrl == 2)
    val = InterpolateQuartic (X, Y, x, index,   index+1, index+2, index+3, order);
  else if (cntrl == 3)
    val = InterpolateQuartic (X, Y, x, index-2, index-1, index,   index+1, order);
  else
    {
      printf ("FLUX::Interpolate: Error - cntrl = %1d\n", cntrl);
      exit (1);
    }
  
  return val;
}

double Flux::InterpolateCubic (double* X, double* Y, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = (x - X[i0]) * (x - X[i2]) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = (x - X[i0]) * (x - X[i1]) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i1]) + (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = ((x - X[i0]) + (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = ((x - X[i0]) + (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else if (order == 2)
    {
      double s0 = 2. /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = 2. /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = 2. /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2];
    }
  else
    {
      printf ("FLUX::InterpolateCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQuartic (double* X, double* Y, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = (x - X[i0]) * (x - X[i2]) * (x - X[i3]) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = (x - X[i0]) * (x - X[i1]) * (x - X[i3]) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = ((x - X[i2]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = ((x - X[i1]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else if (order == 2)
    {
      double s0 = 2. * ((x - X[i1]) + (x - X[i2]) + (x - X[i3])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = 2. * ((x - X[i0]) + (x - X[i2]) + (x - X[i3])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i3])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i2])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * Y[i0] + s1 * Y[i1] + s2 * Y[i2] + s3 * Y[i3];
    }
  else
    {
      printf ("FLUX::InterpolateQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ########################################################
// 1D interpolation function with nonuniform monotonic grid
// order = 0: Y(x)
// order = 1: dY/dx
// order = 2: d^2Y/dx^2
// ########################################################
double Flux::Interpolate1 (int I, double* X, gsl_matrix* Y, int j, double x, int order)
{
  int index, cntrl = 0;
  
  if (x < X[0])
    {
      index = 0;
      cntrl = 2;
    }
  else if (x >= X[I-1])
    {
      index = I - 2;
      cntrl = 3;
    }
  else
    {
      for (int i = 0; i < I-1; i++)
	if (x >= X[i] && x < X[i+1])
	  {
	    index = i;

	    if (index == 0)
	      cntrl = 2;
	    else if (index == I-2)
	      cntrl = 3;
	    else
	      cntrl = 1;

	    break;
	  }
    }

  double val;
  if (cntrl == 1)
    val = InterpolateQuartic1 (X, Y, j, x, index-1, index,   index+1, index+2, order);
  else if (cntrl == 2)
    val = InterpolateQuartic1 (X, Y, j, x, index,   index+1, index+2, index+3, order);
  else if (cntrl == 3)
    val = InterpolateQuartic1 (X, Y, j, x, index-2, index-1, index,   index+1, order);
  else
    {
      printf ("FLUX::Interpolate2: Error - cntrl = %1d\n", cntrl);
      exit (1);
    }
  
  return val;
}

double Flux::InterpolateCubic1 (double* X, gsl_matrix* Y, int j, double x, int i0, int i1, int i2, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = (x - X[i0]) * (x - X[i2]) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = (x - X[i0]) * (x - X[i1]) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
      
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 * gsl_matrix_get (Y, i1, j) + s2 * gsl_matrix_get (Y, i2, j);
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i1]) + (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = ((x - X[i0]) + (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = ((x - X[i0]) + (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 * gsl_matrix_get (Y, i1, j) + s2 * gsl_matrix_get (Y, i2, j);
    }
  else if (order == 2)
    {
      double s0 = 2. /(X[i0] - X[i1]) /(X[i0] - X[i2]);
      double s1 = 2. /(X[i1] - X[i0]) /(X[i1] - X[i2]);
      double s2 = 2. /(X[i2] - X[i0]) /(X[i2] - X[i1]);
  
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 * gsl_matrix_get (Y, i1, j) + s2 * gsl_matrix_get (Y, i2, j);
    }
  else
    {
      printf ("FLUX::InterpolateCubic1: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQuartic1 (double* X, gsl_matrix* Y, int j, double x, int i0, int i1, int i2, int i3, int order)
{  
  double val;

  if (order == 0)
    {
      double s0 = (x - X[i1]) * (x - X[i2]) * (x - X[i3]) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = (x - X[i0]) * (x - X[i2]) * (x - X[i3]) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = (x - X[i0]) * (x - X[i1]) * (x - X[i3]) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = (x - X[i0]) * (x - X[i1]) * (x - X[i2]) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 *  gsl_matrix_get (Y, i1, j) + s2 *  gsl_matrix_get (Y, i2, j) + s3 *  gsl_matrix_get (Y, i3, j);
    }
  else if (order == 1)
    {
      double s0 = ((x - X[i2]) * (x - X[i3]) + (x - X[i1]) * (x - X[i3]) + (x - X[i1]) * (x - X[i2])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = ((x - X[i2]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i2])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = ((x - X[i1]) * (x - X[i3]) + (x - X[i0]) * (x - X[i3]) + (x - X[i0]) * (x - X[i1])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = ((x - X[i1]) * (x - X[i2]) + (x - X[i0]) * (x - X[i2]) + (x - X[i0]) * (x - X[i1])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 *  gsl_matrix_get (Y, i1, j) + s2 *  gsl_matrix_get (Y, i2, j) + s3 *  gsl_matrix_get (Y, i3, j);
    }
  else if (order == 2)
    {
      double s0 = 2. * ((x - X[i1]) + (x - X[i2]) + (x - X[i3])) /(X[i0] - X[i1]) /(X[i0] - X[i2]) /(X[i0] - X[i3]);
      double s1 = 2. * ((x - X[i0]) + (x - X[i2]) + (x - X[i3])) /(X[i1] - X[i0]) /(X[i1] - X[i2]) /(X[i1] - X[i3]);
      double s2 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i3])) /(X[i2] - X[i0]) /(X[i2] - X[i1]) /(X[i2] - X[i3]);
      double s3 = 2. * ((x - X[i0]) + (x - X[i1]) + (x - X[i2])) /(X[i3] - X[i0]) /(X[i3] - X[i1]) /(X[i3] - X[i2]);
      
      val = s0 * gsl_matrix_get (Y, i0, j) + s1 *  gsl_matrix_get (Y, i1, j) + s2 *  gsl_matrix_get (Y, i2, j) + s3 *  gsl_matrix_get (Y, i3, j);
    }
  else
    {
      printf ("FLUX::InterpolateQuartic1: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ###########################################################
// Function to interpolate on uniform 2D grid
// Defaults to quartic interpolation
// Uses cubic/quardatic iterpolation close to grid boundaries
//
// order = 0: Psi
// order = 1: Psi_r
// order = 2: Psi_z
// order = 3; Psi_rr
// order = 4; Psi_rz
// order = 5; Psi_zz
//
// cubic is flag for only cubic interpolation
//
// ###########################################################
double Flux::Interpolate (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX, int order, int cubic)
{
  int i0, ic = 0;
  
  if (RR < Rpts[0])
    {
      i0 = 0;
      ic = 2;
    }
  else if (RR >= Rpts[NRpts-1])
    {
      i0 = NRpts - 2;
      ic = 3;
    }
  else
    {
      for (int i = 0; i < NRpts-1; i++)
	if (RR >= Rpts[i] && RR < Rpts[i+1])
	  {
	    i0 = i;

	    if (i0 == 0)
	      ic = 2;
	    else if (i0 == NRpts-2)
	      ic = 3;
	    else
	      ic = 1;

	    break;
	  }
    }

  int j0, jc = 0;
  
  if (ZZ < Zpts[0])
    {
      j0 = 0;
      jc = 2;
    }
  else if (ZZ >= Zpts[NZpts-1])
    {
      j0 = NZpts - 2;
      jc = 3;
    }
  else
    {
      for (int i = 0; i < NZpts-1; i++)
	if (ZZ >= Zpts[i] && ZZ < Zpts[i+1])
	  {
	    j0 = i;

	    if (j0 == 0)
	      jc = 2;
	    else if (j0 == NZpts-2)
	      jc = 3;
	    else
	      jc = 1;

	    break;
	  }
    }

  double val;
  if (cubic)
    {
      if      (ic == 1 && jc == 1)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0-1, j0,   j0+1, order);
      else if (ic == 1 && jc == 2)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0,   j0+1, j0+2, order);
      else if (ic == 1 && jc == 3)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0-1, j0,   j0+1, order);
      else if (ic == 2 && jc == 1)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0,   i0+1, i0+2, j0-1, j0,   j0+1, order);
      else if (ic == 3 && jc == 1)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0-1, j0,   j0+1, order);
      else if (ic == 2 && jc == 2)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0,   i0+1, i0+2, j0,   j0+1, j0+2, order);
      else if (ic == 2 && jc == 3)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0,   i0+1, i0+2, j0-1, j0,   j0+1, order);
      else if (ic == 3 && jc == 2)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0,   j0+1, j0+2, order);
      else if (ic == 3 && jc == 3)
	val = InterpolateCubicCubic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0,   i0+1, j0-1, j0,   j0+1, order);
       else
	{
	  printf ("FLUX::Interpolate: Error - ic = %1d, jc = %1d\n", ic, jc);
	  exit (1);
	}
    }
  else
    {
      if      (ic == 1 && jc == 1)
	val = InterpolateQuarticQuartic (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1, i0+2, j0-1, j0, j0+1, j0+2, order);
      else if (ic == 1 && jc == 2)
	val = InterpolateQuarticCubic   (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1, i0+2,       j0, j0+1, j0+2, order);
      else if (ic == 1 && jc == 3)
	val = InterpolateQuarticCubic   (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1, i0+2, j0-1, j0, j0+1,       order);
      else if (ic == 2 && jc == 1)
	val = InterpolateCubicQuartic   (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX,       i0, i0+1, i0+2, j0-1, j0, j0+1, j0+2, order);
      else if (ic == 3 && jc == 1)
	val = InterpolateCubicQuartic   (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1,       j0-1, j0, j0+1, j0+2, order);
      else if (ic == 2 && jc == 2)
	val = InterpolateCubicCubic     (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX,       i0, i0+1, i0+2,       j0, j0+1, j0+2, order);
      else if (ic == 2 && jc == 3)
	val = InterpolateCubicCubic     (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX,       i0, i0+1, i0+2, j0-1, j0, j0+1,       order);
      else if (ic == 3 && jc == 2)
	val = InterpolateCubicCubic     (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1,             j0, j0+1, j0+2, order);
      else if (ic == 3 && jc == 3)
	val = InterpolateCubicCubic     (RR, ZZ, NRpts, NZpts, Rpts, Zpts, XX, i0-1, i0, i0+1,       j0-1, j0, j0+1,       order);
      else
	{
	  printf ("FLUX::Interpolate: Error - ic = %1d, jc = %1d\n", ic, jc);
	  exit (1);
	}
    }
  
  return val;
}

double Flux::InterpolateCubicCubic (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX, int i0, int i1, int i2, int j0, int j1, int j2, int order)
{
  double val;

  double dR  = RPTS[1] - RPTS[0];
  double dZ  = ZPTS[1] - ZPTS[0];
  double dR2 = dR*dR;
  double dZ2 = dZ*dZ;

  double RR0 = Rpts[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;

  double ZZ0 = Zpts[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
      
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
       
      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 3)
    {
      double r0 = 2. /(2.*dR2);
      double r1 = 2. /(-  dR2);
      double r2 = 2. /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);
       
      double z0 = 2. /(2.*dZ2);
      double z1 = 2. /(-  dZ2);
      double z2 = 2. /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2;
    }
  else
    {
      printf ("FLUX::InterpolateCubicCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQuarticCubic (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int order)
{
  double val;

  double dR  = Rpts[1] - Rpts[0];
  double dZ  = Zpts[1] - Zpts[0];
  double dR2 = dR*dR;
  double dZ2 = dZ*dZ;
  double dR3 = dR*dR2;
  
  double RR0 = Rpts[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;
  double RR3 = RR2 + dR;

  double ZZ0 = Zpts[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
      
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);
 
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);
      
      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
       
      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);
   
      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
   else if (order == 3)
    {
      double r0 = 2. * ((RR - RR1) + (RR - RR2) + (RR - RR3)) /(-6.*dR3);
      double r1 = 2. * ((RR - RR0) + (RR - RR2) + (RR - RR3)) /(+2.*dR3);
      double r2 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR3)) /(-2.*dR3);
      double r3 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR2)) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) /(2.*dZ2);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) /(-  dZ2);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) /(2.*dZ2);
 
      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = ((ZZ - ZZ1) + (ZZ - ZZ2)) /(2.*dZ2);
      double z1 = ((ZZ - ZZ0) + (ZZ - ZZ2)) /(-  dZ2);
      double z2 = ((ZZ - ZZ0) + (ZZ - ZZ1)) /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val + r3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
       
      double z0 = 2. /(2.*dZ2);
      double z1 = 2. /(-  dZ2);
      double z2 = 2. /(2.*dZ2);

      double val0 = z0 * XX (i0, j0) + z1 * XX (i0, j1) + z2 * XX (i0, j2);
      double val1 = z0 * XX (i1, j0) + z1 * XX (i1, j1) + z2 * XX (i1, j2);
      double val2 = z0 * XX (i2, j0) + z1 * XX (i2, j1) + z2 * XX (i2, j2);
      double val3 = z0 * XX (i3, j0) + z1 * XX (i3, j1) + z2 * XX (i3, j2);

      val = r0 * val0 + r1 * val1 + r2 * val2 + r3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolateQuarticCubic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateCubicQuartic (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX, int i0, int i1, int i2, int j0, int j1, int j2, int j3, int order)
{
  double val;

  double dR  = Rpts[1] - Rpts[0];
  double dZ  = Zpts[1] - Zpts[0];
  double dR2 = dR*dR;
  double dZ2 = dZ*dZ;
  double dZ3 = dZ*dZ2;
  
  double RR0 = Rpts[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;

  double ZZ0 = Zpts[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;
  double ZZ3 = ZZ2 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);
 
      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 3)
    {
      double r0 = 2. /(2.*dR2);
      double r1 = 2. /(-  dR2);
      double r2 = 2. /(2.*dR2);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR1) + (RR - RR2)) /(2.*dR2);
      double r1 = ((RR - RR0) + (RR - RR2)) /(-  dR2);
      double r2 = ((RR - RR0) + (RR - RR1)) /(2.*dR2);
      
      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) /(2.*dR2);
      double r1 = (RR - RR0) * (RR - RR2) /(-  dR2);
      double r2 = (RR - RR0) * (RR - RR1) /(2.*dR2);

      double z0 = 2. * ((ZZ - ZZ1) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(-6.*dZ3);
      double z1 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(+2.*dZ3);
      double z2 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ3)) /(-2.*dZ3);
      double z3 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ2)) /(+6.*dZ3);
 
      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolateCubicQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

double Flux::InterpolateQuarticQuartic (double RR, double ZZ, int NRpts, int NZpts, double* Rpts, double* Zpts, Array<double,2> XX, int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, int order)
{
  double val;

  double dR  = Rpts[1] - Rpts[0];
  double dZ  = Zpts[1] - Zpts[0];
  double dR2 = dR*dR;
  double dZ2 = dZ*dZ;
  double dR3 = dR*dR2;
  double dZ3 = dZ*dZ2;
  
  double RR0 = Rpts[i0];
  double RR1 = RR0 + dR;
  double RR2 = RR1 + dR;
  double RR3 = RR2 + dR;

  double ZZ0 = Zpts[j0];
  double ZZ1 = ZZ0 + dZ;
  double ZZ2 = ZZ1 + dZ;
  double ZZ3 = ZZ2 + dZ;

  if (order == 0)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 1)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 2)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);
      
      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 3)
    {
      double r0 = 2. * ((RR - RR1) + (RR - RR2) + (RR - RR3)) /(-6.*dR3);
      double r1 = 2. * ((RR - RR0) + (RR - RR2) + (RR - RR3)) /(+2.*dR3);
      double r2 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR3)) /(-2.*dR3);
      double r3 = 2. * ((RR - RR0) + (RR - RR1) + (RR - RR2)) /(+6.*dR3);
  
      double z0 = (ZZ - ZZ1) * (ZZ - ZZ2) * (ZZ - ZZ3) /(-6.*dZ3);
      double z1 = (ZZ - ZZ0) * (ZZ - ZZ2) * (ZZ - ZZ3) /(+2.*dZ3);
      double z2 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ3) /(-2.*dZ3);
      double z3 = (ZZ - ZZ0) * (ZZ - ZZ1) * (ZZ - ZZ2) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 4)
    {
      double r0 = ((RR - RR2) * (RR - RR3) + (RR - RR1) * (RR - RR3) + (RR - RR1) * (RR - RR2)) /(-6.*dR3);
      double r1 = ((RR - RR2) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR2)) /(+2.*dR3);
      double r2 = ((RR - RR1) * (RR - RR3) + (RR - RR0) * (RR - RR3) + (RR - RR0) * (RR - RR1)) /(-2.*dR3);
      double r3 = ((RR - RR1) * (RR - RR2) + (RR - RR0) * (RR - RR2) + (RR - RR0) * (RR - RR1)) /(+6.*dR3);

      double z0 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ1) * (ZZ - ZZ2)) /(-6.*dZ3);
      double z1 = ((ZZ - ZZ2) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ2)) /(+2.*dZ3);
      double z2 = ((ZZ - ZZ1) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ3) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(-2.*dZ3);
      double z3 = ((ZZ - ZZ1) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ2) + (ZZ - ZZ0) * (ZZ - ZZ1)) /(+6.*dZ3);

      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else if (order == 5)
    {
      double r0 = (RR - RR1) * (RR - RR2) * (RR - RR3) /(-6.*dR3);
      double r1 = (RR - RR0) * (RR - RR2) * (RR - RR3) /(+2.*dR3);
      double r2 = (RR - RR0) * (RR - RR1) * (RR - RR3) /(-2.*dR3);
      double r3 = (RR - RR0) * (RR - RR1) * (RR - RR2) /(+6.*dR3);

      double z0 = 2. * ((ZZ - ZZ1) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(-6.*dZ3);
      double z1 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ2) + (ZZ - ZZ3)) /(+2.*dZ3);
      double z2 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ3)) /(-2.*dZ3);
      double z3 = 2. * ((ZZ - ZZ0) + (ZZ - ZZ1) + (ZZ - ZZ2)) /(+6.*dZ3);
 
      double val0 = r0 * XX (i0, j0) + r1 * XX (i1, j0) + r2 * XX (i2, j0) + r3 * XX (i3, j0);
      double val1 = r0 * XX (i0, j1) + r1 * XX (i1, j1) + r2 * XX (i2, j1) + r3 * XX (i3, j1);
      double val2 = r0 * XX (i0, j2) + r1 * XX (i1, j2) + r2 * XX (i2, j2) + r3 * XX (i3, j2);
      double val3 = r0 * XX (i0, j3) + r1 * XX (i1, j3) + r2 * XX (i2, j3) + r3 * XX (i3, j3);

      val = z0 * val0 + z1 * val1 + z2 * val2 + z3 * val3;
    }
  else
    {
      printf ("FLUX::InterpolateQuarticQuartic: Error - order = %1d\n", order);
      exit (1);
    }

  return val;
}

// ###############################
// Function to evaluate Psi (R, Z)
// ###############################
double Flux::GetPsi (double r, double z)
{
  return Interpolate (r, z, NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 0, 0);
}

// ###################################
// Function to evaluate dPsi/dR (R, Z)
// ###################################
double Flux::GetPsiR (double r, double z)
{
  return Interpolate (r, z, NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 1, 0);
}

// ###################################
// Function to evaluate dPsi/dZ (R, Z)
// ###################################
double Flux::GetPsiZ (double r, double z)
{
  return Interpolate (r, z, NRPTS, NZPTS, RPTS, ZPTS, PSIARRAY, 2, 0);
}


