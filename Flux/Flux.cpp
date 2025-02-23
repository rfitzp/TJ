// Flux.cpp

#include "Flux.h"

// ###########
// Constructor
// ###########
Flux::Flux ()
{
   // -------------------------------------------
  // Ensure that directory ../Outputs/Flux exits
  // -------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Flux"))
    {
      exit (1);
    }
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
  // --------------------------------------
  // Set default values of input parameters
  // --------------------------------------
  NPSI   = 256;
  PACK   = 1.;
  NTHETA = 513;
  H0     = 1.e-6;
  ACC    = 1.e-14;

  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/Flux.json";
  json   JSONData     = ReadJSONFile (JSONFilename);
  
  NPSI   = JSONData["NPSI"]  .get<int>   ();
  PACK   = JSONData["PACK"]  .get<double>();
  NTHETA = JSONData["NTHETA"].get<int>   ();
  H0     = JSONData["H0"]    .get<double>();
  ACC    = JSONData["ACC"]   .get<double>();

  // ------------
  // Sanity check
  // ------------
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

  // -----------------------------
  // Output calculation parameters
  // -----------------------------
  printf ("\n");
  printf ("Class FLUX::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Input Parameters (from Inputs/Flux.json):\n");
  printf ("NPSI = %4d        NTHETA = %4d       PACK = %10.3e\n",
	  NPSI, NTHETA, PACK);
  printf ("H0   = %10.3e  ACC    = %10.3e\n",
	  H0, ACC);
 }

// ###################################################
// Function to input gFile data and output Stage1 data
// ###################################################
void Flux::Stage1 ()
{
  // Read gFile
  gFileRead ();
}

