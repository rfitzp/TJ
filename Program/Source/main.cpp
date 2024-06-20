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
  Equilibrium equilibrium;
  equilibrium.Solve ();
  
  // ...........................................................
  // Call class TJ to calculate tearing stability of equilibrium
  // ...........................................................
  TJ tj;
  tj.Solve ();

  return 0;
}

