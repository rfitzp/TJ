// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program StartUp
// See StartUp.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "StartUp.h"

int main (int argc, char* argv[])
{
  printf ("---------------\nProgram StartUp\n---------------\n");
  clock_t begin = clock ();

  StartUp startup;
  startup.Solve ();
  
  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);
}

