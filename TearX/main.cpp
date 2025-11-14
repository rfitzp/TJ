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

  printf ("\nOption ? (Single calculation = 1 q0 scan = 2 nu scan = 3) ");
  scanf ("%d", &calc);
  if (calc < 1 || calc > 3)
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
      tearx.Scanq0 ();
    }
  else if (calc == 3)
    {
      tearx.Scannu ();
    }
    
  return 0;
}
