// Flux.cpp

// PROGRAM ORGANIZATION:
//
//       Flux:: Flux          ()
//       Flux:: ~Flux         ()       
// void  Flux:: Solve         ()
// void  Flux:: SetParameters ()
// FILE* Flux:: OpenFilew     (char* filename)
// FILE* Flux:: OpenFiler     (char* filename)
// FILE* Flux:: OpenFilea     (char* filename)

#include "Flux.h"

// ###########
// Constructor
// ###########
Flux::Flux ()
{
}

// ##########
// Destructor
// ##########
Flux::~Flux ()
{
}

// #########################
// Function to solve problem
// #########################
void Flux::Solve ()
{
  // Set global parameters
  SetParameters ();

  // Input gFile data and output Stage1 data.
  // Stage1 data output to directory Outputs/Stage1.
  Stage1 ();
  
  // Input Stage1 data and output Stage2 data.
  // Stage2 data output to directory Outputs/Stage2.
  Stage2 ();
}

// #################################
// Function to set global parameters
// #################################
void Flux::SetParameters ()
{
  // Set default values of input parameters
  NPSI    = 256;
  PACK    = 1.;
  NTHETA  = 512;
  H0      = 1.e-6;
  ACC     = 1.e-14;

  // Read namelist file Inputs/Flux.nml
  NameListRead (&NPSI, &PACK, &NTHETA, &H0, &ACC);

  // Sanity check
  if (NPSI < 1)
    {
      printf ("FLUX::SetParameters: Error - NPSI must be positive\n");
      exit (1);
    }
  if (NTHETA < 1)
    {
      printf ("FLUX::SetParameters: Error NTHETA must be positive\n");
      exit (1);
    }
  if (NTHETA%2 == 0)
    {
      printf ("FLUX::SetParameters: Error - NTHETA must be odd\n");
      exit (1);
    }
  if (H0 <= 0.)
    {
      printf ("FLUX::SetParameters: Error - H0 must be positive\n");
      exit (1);
    }
  if (ACC <= 0.)
    {
      printf ("FLUX::SetParameters: Error - ACC must be positive\n");
      exit (1);
    }
  
   // Output calculation parameters
   printf ("Input Parameters (from Inputs/Flux.nml):\n");
   printf ("NPSI = %4d  NTHETA = %4d  PACK = %10.3e \n",
	   NPSI, NTHETA, PACK);
    printf ("H0   = %10.3e  ACC = %10.3e\n",
	   H0, ACC);
 }

// #####################################
// Function to open new file for writing
// #####################################
FILE* Flux::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* Flux::OpenFiler (char* filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ############################################
// Function to open existing file for appending
// ############################################
FILE* Flux::OpenFilea (char* filename)
{
  FILE* file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("FLUX::OpenFilea: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// #################################
// Function to call operating system
// #################################
void Flux::CallSystem (char* command)
{
  if (system (command) != 0)
    {
      printf ("FLUX: Operating system call error executing %s\n", command);
      exit (1);
    }
}
