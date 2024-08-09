// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TJ
// See Equilibium.h, TJ.h, and Layer.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Equilibrium.h"
#include "TJ.h"
#include "Layer.h"

int main ()
{
  printf ("----------\nProgram TJ\n----------\n");
  clock_t begin = clock ();

  // .............................................................................
  // Call class Equilibrium to construct aspect-ratio expanded tokamak equilibrium
  // .............................................................................
  {
    Equilibrium equilibrium;
    equilibrium.Solve ();
  }
  
  // ...........................................................
  // Call class TJ to calculate tearing stability of equilibrium
  // ...........................................................
  {
    TJ tj;
    tj.Solve ();
  }

  // ...............................................................
  // Call class Layer to calculate growth-rates and real frequencies
  // ...............................................................
  {
    Layer layer;
    layer.Solve ();
  }

  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);

  return 0;
}

