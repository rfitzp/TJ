// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Layer
// See Layer.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Layer.h"

// Pointer to target function for passing to GSL multidimensional root finding algorithm
int pTarget (const gsl_vector* x, void* params, gsl_vector* f)
{
  Layer layer = *(Layer*) params;

  const double x1 = gsl_vector_get (x, 0);
  const double x2 = gsl_vector_get (x, 1);

  double F1, F2;

  layer.Target (x1, x2, F1, F2);

  gsl_vector_set (f, 0, F1);
  gsl_vector_set (f, 1, F2);

  return GSL_SUCCESS;
}

int main (int argc, char** argv)
{
  // ..................
  // Call program Layer
  // ..................
  Layer layer;
  layer.Solve (1);

  return 0;
}
