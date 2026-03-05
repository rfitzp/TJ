// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program NonLinear
// See NonLinear.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "NonLinear.h"

int main ()
{
  // ..............
  // Select machine
  // ..............
  printf ("\nMachine ? (1 = RFP, 2 = Tokamak) ");
  int _machine;
  scanf ("%d", &_machine);
  if (_machine < 1 || _machine > 2)
    {
      printf ("Error - Unknown machine\n");
      exit (1);
    }
    
  // ................
  // Calculation menu
  // ................
  int calc;
  printf ("\nOption ? (Single calculation = 1, beta0 scan  = 2, W scan = 3, m=0 = 4, m=1 = 5, m=1/m=0 6) ");
  scanf ("%d", &calc);
  if (calc < 1 || calc > 6)
    {
      printf ("Error - Unknown option\n");
      exit (1);
    }

  if (calc == 1)
    {
      NonLinear nonlinear (_machine);
      double tau, taup, Delta1, Delta2, Delta3;
      nonlinear.Solve (1, tau, taup, Delta1, Delta2, Delta3);
      printf ("tau (no pressure jumps) = %11.4e  tau (with pressure jumps) = %11.4e\n", tau, taup);
     }
  else if (calc == 2)
    {
      NonLinear nonlinear (_machine);
      double b0s, b0e;
      int    Nb;
      printf ("\nbeta0_start, beta0_end, Nbeta0 ??? ");
      scanf ("%lf %lf %d", &b0s, &b0e, &Nb);
      if (b0s >= 0. && b0e > 0. && Nb > 2)
	{
	  printf ("\n");
	  FILE* file = fopen ("Output/Bscan.out", "w");
	  for (int i = 0; i <= Nb; i++)
	    {
	      double _beta0 = b0s + (b0e - b0s ) * double (i) /double (Nb);
	      double tau, taup, Delta1, Delta2, Delta3;
	      nonlinear.Setbeta0 (_beta0);
	      nonlinear.Solve (0, tau, taup, Delta1, Delta2, Delta3);

	      printf ("beta0 = %11.4e  tau = %11.4e  taup = %11.4e\n", _beta0, tau, taup);

	      if (_machine == 1)
		{
		  if (taup > 0.)
		    fprintf (file, "%e %e %e\n", _beta0, log10(tau), log10(taup));
		  else if (_beta0 < 0.04)
		    fprintf (file, "%e %e %e\n", _beta0, log10(tau), log10(-taup));
		}
	      else
		fprintf (file, "%e %e %e\n", _beta0, tau, taup);
	    }
	  fclose (file);
	}
      else
	{
	  printf ("Bad scan parameters\n");
	  exit (1);
	}
    }
  else if (calc == 3)
    {
      NonLinear nonlinear (_machine);
      double Ws, We;
      int    NW;
      printf ("\nW_start, W_end, NW ??? ");
      scanf ("%lf %lf %d", &Ws, &We, &NW);
      if (Ws >= 0. && We > 0. && NW > 2)
	{
	  printf ("\n");
	  FILE* file = fopen ("Output/Wscan.out", "w");
	  for (int i = 0; i <= NW; i++)
	    {
	      double _W = Ws + (We - Ws ) * double (i) /double (NW);
	      double tau, taup, Delta1, Delta2, Delta3;
	      nonlinear.SetW (_W);
	      nonlinear.Solve (0, tau, taup, Delta1, Delta2, Delta3);

	      printf ("W = %11.4e  tau = %11.4e  taup = %11.4e\n", _W, tau, taup);

	      fprintf (file, "%e %e %e\n", _W, tau, taup);
	    }
	  fclose (file);
	}
      else
	{
	  printf ("Bad scan parameters\n");
	  exit (1);
	}
    }
  else if (calc == 4)
    {
      NonLinear nonlinear (_machine);
      int _mpol = 0;
      double _beta0;
      printf ("beta0 ?? ");
      scanf ("%lf", &_beta0);
      nonlinear.Setbeta0 (_beta0);
      FILE* file = fopen ("Output/m0.out", "w");
      for (int _ntor = 1; _ntor < 100; _ntor++)
	{
	  double Delta, Psirv;
	  nonlinear.GetDelta (_mpol, _ntor, Delta, Psirv);
	  printf ("mpol = %3d  ntor = %3d  Delta = %11.4e  Psi(r_v) = %11.4e\n", _mpol, _ntor, Delta, Psirv);
	  if (Delta > 0.)
	    fprintf (file, "%11.4e %11.4e %11.4e\n", double (_ntor), Delta, Psirv);
	  else
	    {
	      fprintf (file, "%11.4e %11.4e %11.4e\n", double (_ntor), Delta, Psirv);
	      break;
	    }
	}
      fclose (file);
    }
  else if (calc == 5)
    {
      NonLinear nonlinear (_machine);
      int _mpol = 1;
      double _beta0;
      printf ("beta0 ?? ");
      scanf ("%lf", &_beta0);
      nonlinear.Setbeta0 (_beta0);
      FILE* file = fopen ("Output/m1.out", "w");
      for (int _ntor = 8; _ntor < 100; _ntor++)
	{
	  double Delta, Psirv;
	  nonlinear.GetDelta (_mpol, _ntor, Delta, Psirv);
	  printf ("mpol = %3d  ntor = %3d  Delta = %11.4e  Psi(r_v) = %11.4e\n", _mpol, _ntor, Delta, Psirv);
	  if (Delta > 0.)
	    fprintf (file, "%11.4e %11.4e %11.4e\n", double (_ntor), Delta, Psirv);
	  else
	    {
	      fprintf (file, "%11.4e %11.4e %11.4e\n", double (_ntor), Delta, Psirv);
	      break;
	    }
	}
      fclose (file);
    }
  else
    {
      NonLinear nonlinear (_machine);
      double _beta0;
      printf ("beta0 ?? ");
      scanf ("%lf", &_beta0);
      nonlinear.Setbeta0 (_beta0);
      double tau, taup, Delta1, Delta2, Delta3;
      FILE* file = fopen ("Output/m10.out", "w");
      for (int _ntor = 8; _ntor < 100; _ntor++)
	{
	  for (int _k = 1; _k < 100; _k++)
	    {
	      nonlinear.SetMode (_ntor, _k);
	      nonlinear.Solve (0, tau, taup, Delta1, Delta2, Delta3);

	      printf ("(n, k) = (%3d, %3d)  Delta = (%11.4e, %11.4e, %11.4e)  J = %11.4e\n",
		      _ntor, _k, Delta1, Delta2, Delta3, taup);

	      if (Delta1 > 0. && Delta2 > 0. && Delta3 > 0.)
		fprintf (file, "%3d %3d %11.4e %11.4e %11.4e %11.4e\n",
			 _ntor, _k, Delta1, Delta2, Delta3, taup);
	      
	      if (Delta3 < 0.) break;
	    }
	  if (Delta1 < 0.) break;
	}
      fclose (file);
    }

  return 0;
}
