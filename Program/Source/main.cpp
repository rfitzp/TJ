// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TJ
// See Equilibium.h and TJ.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Equilibrium.h"
#include "TJ.h"

int main ()
{
  // .............................................................................
  // Call class Equilibrium to construct aspect-ratio expanded tokamak equilibrium
  // .............................................................................
  printf ("----------\nProgram TJ\n----------\n");
  Equilibrium equilibrium;
  clock_t     begin = clock ();
  equilibrium.Solve ();
  
  // ...........................................................
  // Call class TJ to calculate tearing stability of equilibrium
  // ...........................................................
  TJ tj;
  tj.Solve ();

  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);

  return 0;
}

