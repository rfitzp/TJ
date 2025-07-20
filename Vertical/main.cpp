// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Vertical
// See Equilibium.h and Vertical.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Vertical.h"

int main (int argc, char* argv[])
{
  printf ("----------------\nProgram Vertical\n----------------\n");
  clock_t begin = clock ();

  Vertical vertical;
  vertical.Solve ();
  
  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);
}

