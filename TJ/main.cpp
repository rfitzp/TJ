// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TJ
// See Equilibium.h, TJ.h, and Layer.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "TJ.h"

#include "../Layer/Layer.h"

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

int main (int argc, char* argv[])
{
  printf ("----------\nProgram TJ\n----------\n");
  clock_t begin = clock ();

  TJ tj;
  tj.Solve ();
  
  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);
}

