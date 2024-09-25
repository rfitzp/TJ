// main.cpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main function for program TJ
// See Equilibium.h, TJ.h, and Layer.h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Equilibrium.h"
#include "TJ.h"
#include "Layer.h"

void RunProgram (int control, double value);
void read_input (char[], char[]);
void read_data  (char[], char[], int, double[]);

int main (int argc, char* argv[])
{
  // Calling program without argument triggers non-iteractive mode
  if (argc == 1)
    {
      RunProgram (0, 0.);

      exit (0);
    }

  // Calling program with argument triggers interactive mode
  char   in_string[5];
  char*  cursor = (char*) "TJ> ";
  double data[5];

  // Main control loop
  printf ("\n");
  while (1)
    {
      read_input (cursor, in_string);

      if ((in_string[0] == 'q') && (in_string[1] == 'u'))
	exit (0);
      else if ((in_string[0] == 'p') && (in_string[1] == 'c'))
	{
	  double pc_start, pc_end;
	  int    Np;
	  printf ("pc_start pc_end Np ?? ");
	  scanf ("%lf %lf %d", &pc_start, &pc_end, &Np);
	  if (pc_start < 0. || pc_end < 0. || Np < 2)
	    {
	      printf ("Error: Invalid range\n");
	      exit (1);
	    }

	  double pc;
	  for (int i = 0; i <= Np; i++)
	    {
	      pc = pc_start + (pc_end - pc_start) * double (i) /double (Np);

	      RunProgram (1, pc);

	      char command1[100], command2[100], command3[100];
	      sprintf (command1, "mv Plots/Equilibrium.nc Runs/Equilibrium.%04d.nc\n", i);
	      sprintf (command2, "mv Plots/TJ.nc Runs/TJ.%04d.nc\n", i);
	      sprintf (command3, "mv Plots/Layer.nc Runs/Layer.%04d.nc\n", i);

	      system (command1);
	      system (command2);
	      system (command3);
	    }
	}
      else if ((in_string[0] == 'h') && (in_string[1] == 'e'))
	{
	  printf ("\n");
	  printf ("pc ... scan pc\n");
	  printf ("he ... print this message\n");
	  printf ("qu ... exit program\n\n");
	}
    }
}

void RunProgram (int control, double value)
  {
  printf ("----------\nProgram TJ\n----------\n");
  clock_t begin = clock ();

  // .............................................................................
  // Call class Equilibrium to construct aspect-ratio expanded tokamak equilibrium
  // .............................................................................
  {
    Equilibrium equilibrium;

    if (control == 1)
      equilibrium.Setpc (value);

    equilibrium.Setnu ();
    
    equilibrium.Solve ();
  }
  
  // ...........................................................
  // Call class TJ to calculate tearing stability of equilibrium
  // ...........................................................
  {
    TJ tj;
    tj.Solve ();
  }

  // ...............................................................
  // Call class Layer to calculate growth-rates and real frequencies
  // ...............................................................
  {
    Layer layer;
    layer.Solve ();
  }

  clock_t end       = clock ();
  double time_spent = double (end - begin) /double(CLOCKS_PER_SEC);

  printf ("\nNormal termination: Cpu time = %10.3e s\n", time_spent);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function to read line of input from terminal and store it in array input[] of length 5
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void read_input (char cursor[], char input[])
{
  int count = 0; char c;

  printf ("%-s", cursor);
  while ((c = getchar ()) && (c != '\n'))
    {
      input[count] = c;
      if (count < 4) ++count;
    }
  input[count] = '\0';
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function to read n data items from terminal and store them in double array x[]
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void read_data (char cursor[], char message[], int n, double x[])
{
  int count = 1; char c;

  printf ("%-s%-s", cursor, message);
  for (count = 1; count <= n; count++)
    scanf ("%lf", &x[count - 1]);
  while ((c = getchar ()) && (c != '\n')) {}
}
