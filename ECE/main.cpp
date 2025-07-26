// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program ECE
// See ECE.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "ECE.h"

int main (int argc, char* argv[])
{
  printf ("-----------\nProgram ECE\n-----------\n");
  clock_t begin = clock ();

  ECE ece;
  ece.Solve ();
  
  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);
}

