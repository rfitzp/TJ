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
      else if ((in_string[0] == 'e') && (in_string[1] == 'p'))
	{
	  double epsa_start, epsa_end, d_epsa;
	  printf ("epsa_start epsa_end d_epsa ?? ");
	  scanf ("%lf %lf %lf", &epsa_start, &epsa_end, &d_epsa);
	  if (epsa_start <= 0. || epsa_end <= 0.)
	    {
	      printf ("Error: Invalid range\n");
	      exit (1);
	    }

	  double epsa = epsa_start;
	  int    i    = 0;
	  do
	    {
	      RunProgram (3, epsa);

	      char command1[100], command2[100], command3[100];
	      sprintf (command1, "cp ../Outputs/Equilibrium/Equilibrium.nc ../Runs/Equilibrium.%04d.nc\n", i);
	      sprintf (command2, "cp ../Outputs/TJ/TJ.nc ../Runs/TJ.%04d.nc\n", i);
	      //	      sprintf (command3, "mv Plots/Layer.nc Runs/Layer.%04d.nc\n", i);

	      system (command1);
	      system (command2);
	      //	      system (command3);

	      epsa += d_epsa;
	      i++;
	    }
	  while (epsa < epsa_end);
	}
      else if ((in_string[0] == 'p') && (in_string[1] == 'c'))
	{
	  double pc_start, pc_end, d_pc;
	  printf ("pc_start pc_end d_pc ?? ");
	  scanf ("%lf %lf %lf", &pc_start, &pc_end, &d_pc);
	  if (pc_start < 0. || pc_end < 0.)
	    {
	      printf ("Error: Invalid range\n");
	      exit (1);
	    }

	  double pc = pc_start;
	  int    i  = 0;
	  do
	    {
	      RunProgram (4, pc);

	      char command1[100], command2[100], command3[100];
	      sprintf (command1, "cp ../Outputs/Equilibrium/Equilibrium.nc ../Runs/Equilibrium.%04d.nc\n", i);
	      sprintf (command2, "cp ../Outputs/TJ/TJ.nc ../Runs/TJ.%04d.nc\n", i);
	      //	      sprintf (command3, "mv Plots/Layer.nc Runs/Layer.%04d.nc\n", i);

	      system (command1);
	      system (command2);
	      //	      system (command3);

	      pc += d_pc;
	      i++;
	    }
	   while (pc < pc_end);
	}
      else if ((in_string[0] == 'h') && (in_string[1] == '2'))
	{
	  double h2_start, h2_end, d_h2;
	  printf ("h2_start h2_end d_h2 ?? ");
	  scanf ("%lf %lf %lf", &h2_start, &h2_end, &d_h2);
	  if (h2_start < 0. || h2_end < 0.)
	    {
	      printf ("Error: Invalid range\n");
	      exit (1);
	    }

	  double h2 = h2_start;
	  int    i  = 0;
	  do
	    {
	      RunProgram (5, h2);

	      char command1[100], command2[100], command3[100];
	      sprintf (command1, "cp ../Outputs/Equilibrium/Equilibrium.nc ../Runs/Equilibrium.%04d.nc\n", i);
	      sprintf (command2, "cp ../Outputs/TJ/TJ.nc ../Runs/TJ.%04d.nc\n", i);
	      //	      sprintf (command3, "mv Plots/Layer.nc Runs/Layer.%04d.nc\n", i);

	      system (command1);
	      system (command2);
	      //	      system (command3);

	      h2 += d_h2;
	      i++;
	    }
	  while (h2 < h2_end);
	}
      else if ((in_string[0] == 'h') && (in_string[1] == 'e'))
	{
	  printf ("\n");
	  printf ("ep ... scan epsa\n");
	  printf ("h2 ... scan H2\n");
	  printf ("h3 ... scan H3\n");
	  printf ("h4 ... scan H4\n");
	  printf ("pc ... scan pc\n");
	  printf ("qc ... scan qc\n");
	  printf ("qa ... scan qa\n");
	  printf ("v2 ... scan V2\n");
	  printf ("v3 ... scan V3\n");
	  printf ("v4 ... scan V4\n");
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
      equilibrium.Setqc (value);
    else if (control == 2)
      equilibrium.Setqa (value);
    else if (control == 3)
      equilibrium.Setepsa (value);
    else if (control == 4)
      equilibrium.Setpc (value);
    else if (control == 5)
      equilibrium.SetH2 (value);
    else if (control == 6)
      equilibrium.SetV2 (value);
    else if (control == 7)
      equilibrium.SetH3 (value);
    else if (control == 8)
      equilibrium.SetV3 (value);
    else if (control == 9)
      equilibrium.SetH4 (value);
    else if (control == 10)
      equilibrium.SetV4 (value);

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
    layer.Solve (0);
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
