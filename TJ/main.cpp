// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TJ
// See Equilibium.h, TJ.h, and Layer.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "TJ.h"

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

