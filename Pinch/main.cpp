// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program Pinch
// See Pinch.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Pinch.h"

int main ()
{
  Pinch pinch;

  pinch.Solve ();
  
  // ................
  // Calculation menu
  // ................
  /*
  int calc;
  printf ("\nOption ? (Single calculation = 1 q scan = 2 Beta scan = 3 nusscan 4) ");
  scanf ("%d", &calc);
  if (calc < 1 || calc > 4)
    {
      printf ("Error - Unknown option\n");
      exit (1);
    }
  */

  /*
  if (calc == 1)
    {
      Pinch pinch;
      double Theta, F;
      pinch.Solve (Theta, F);
      printf ("Theta = %11.4e  F = %11.4e\n", Theta, F);
    }
  else if (calc == 2)
    {
      Pinch pinch;
      double q0, Theta, F;

      int    N    = 100;
      double qmin = 0.12;
      double qmax = 0.5;
      FILE*  file = fopen ("Output/q0scan000.out", "w");
      pinch.Setbeta0 (0.);
      for (int i = 0; i <= N; i++)
	{
	  q0 = qmin + (qmax - qmin) * double (i) /double (N);
	  pinch.Setq0 (q0);
	  pinch.Solve (Theta, F);
	  printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	  fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	}
      fclose (file);

      file = fopen ("Output/q0scan050.out", "w");
      pinch.Setbeta0 (0.05);
      for (int i = 0; i <= N; i++)
	{
	  q0 = qmin + (qmax - qmin) * double (i) /double (N);
	  pinch.Setq0 (q0);
	  pinch.Solve (Theta, F);
	  printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	  fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	}
      fclose (file);

      file = fopen ("Output/q0scan100.out", "w");
      pinch.Setbeta0 (0.1);
      for (int i = 0; i <= N; i++)
	{
	  q0 = qmin + (qmax - qmin) * double (i) /double (N);
	  pinch.Setq0 (q0);
	  pinch.Solve (Theta, F);
	  printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	  fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	}
      fclose (file);
    }
  else if (calc == 3)
    {
      Pinch pinch;
      double beta0, Theta, F;

      int    N       = 100;
      double betamax = 0.1;
      for (int i = 0; i <= N; i++)
	{
	  beta0 = betamax * double (i) /double (N);
	  pinch.Setbeta0 (beta0);
	  pinch.Solve (Theta, F);
	  printf ("beta0 = %11.4e  Theta = %11.4e  F = %11.4e\n", beta0, Theta, F);
	}
    }
   else if (calc == 4)
     {
       Pinch pinch;
       double q0, Theta, F;

       int    N    = 100;
       double qmin = 0.125;
       double qmax = 0.5;
       FILE*  file = fopen ("Output/nusscan100.out", "w");
       pinch.Setbeta0 (0.05);
       pinch.Setnus (1.);
       for (int i = 0; i <= N; i++)
	 {
	   q0 = qmin + (qmax - qmin) * double (i) /double (N);
	   pinch.Setq0 (q0);
	   pinch.Solve (Theta, F);
	   printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	   fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	 }
       fclose (file);
       
       file = fopen ("Output/nusscan150.out", "w");
       pinch.Setnus (1.5);
       for (int i = 0; i <= N; i++)
	 {
	   q0 = qmin + (qmax - qmin) * double (i) /double (N);
	   pinch.Setq0 (q0);
	   pinch.Solve (Theta, F);
	   printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	   fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	 }
       fclose (file);
       
       file = fopen ("Output/nusscan200.out", "w");
       pinch.Setnus (2.0);
       for (int i = 0; i <= N; i++)
	 {
	   q0 = qmin + (qmax - qmin) * double (i) /double (N);
	   pinch.Setq0 (q0);
	   pinch.Solve (Theta, F);
	   printf ("q0 = %11.4e  Theta = %11.4e  F = %11.4e\n", q0, Theta, F);
	   fprintf (file, "%11.4e %11.4e %11.4e\n", q0, Theta, F);
	 }
       fclose (file);
     }
  */
  
  return 0;
}
