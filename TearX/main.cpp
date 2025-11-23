// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TearX
// See TearX.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "TearX.h"

int main (int argc, char** argv)
{
  // ................
  // Calculation menu
  // ................
  int   calc;
  TearX tearx;

  printf ("\nOption ? (Single calculation = 1 nu scan = 2) ");
  scanf ("%d", &calc);
  if (calc < 1 || calc > 2)
    {
      printf ("Error - Unknown option\n");
      exit (1);
    }

  if (calc == 1)
    {
      tearx.Solve (0);
    }
  else if (calc == 2)
    {
      tearx.Scannu ();
    }
    
  return 0;
}
