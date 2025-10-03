// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Pinch
// See Pinch.h
//
// Command line options:
//
// 1 - Scan toroidal mode number
// 2 - Scan poloidal mode number
// 3 - Scan wall position
// 4 - Scan central beta
// 5 - Scan q0
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Pinch.h"

int main (int argc, char* argv[])
{
  Pinch pinch;

  if (argc == 1)
    pinch.Solve (1, -1.);
  else
    {
      int n = atoi (argv[1]);
      if (n < 1 || n > 5)
	{
	  printf ("Unknown command line option\n");
	  exit (1);
	}
      else
	pinch.Scan (n);
    }
   
  return 0;
}
