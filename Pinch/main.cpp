// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Pinch
// See Pinch.h
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
      if (n < 1|| n > 4)
	{
	  printf ("Unknown command line option\n");
	  exit (1);
	}
      else
	pinch.Scan (n);
    }
   
  return 0;
}
