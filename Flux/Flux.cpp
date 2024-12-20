// Flux.cpp

// PROGRAM ORGANIZATION:
//
//       Flux:: Flux          ()
//       Flux:: ~Flux         ()       
// void  Flux:: Solve         ()
// void  Flux:: SetParameters ()
// void  Flux:: Stage1        ()
// json  Flux:: ReadJSONFile  (const string& filename)
// FILE* Flux:: OpenFilew     (char* filename)
// FILE* Flux:: OpenFiler     (char* filename)
// FILE* Flux:: OpenFilea     (char* filename)

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

// ########################################
// Function to strip comments from a string
// ########################################
string Flux::stripComments (const string& input)
{
  stringstream result;
  bool         inSingleLineComment = false;
  bool         inMultiLineComment  = false;

  for (size_t i = 0; i < input.size(); ++i)
    {
      // Start of single-line comment (//)
      if (!inMultiLineComment && input[i] == '/' && i + 1 < input.size() && input[i + 1] == '/')
	{
	  inSingleLineComment = true;
	  i++; 
	}
      // Start of multi-line comment (/* ... */)
      else if (!inSingleLineComment && input[i] == '/' && i + 1 < input.size() && input[i + 1] == '*')
	{
	  inMultiLineComment = true;
	  i++; 
	}
      // End of single-line comment
      else if (inSingleLineComment && input[i] == '\n')
	{
	  inSingleLineComment = false;
	  result << input[i];
	}
      // End of multi-line comment
      else if (inMultiLineComment && input[i] == '*' && i + 1 < input.size() && input[i + 1] == '/')
	{
	  inMultiLineComment = false;
	  i++; 
	}
      // Regular characters outside comments
      else if (!inSingleLineComment && !inMultiLineComment)
	{
	  result << input[i];
	}
    }
  
  return result.str();
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

// ##########################
// Function to read JSON file
// ##########################
json Flux::ReadJSONFile (const string& filename)
{
  ifstream JSONFile (filename);
  json     JSONData;

  if (JSONFile.is_open ())
    {
      try
	{
	  // Strip any comments from JSON file
	  stringstream buffer;
	  buffer << JSONFile.rdbuf ();
	  JSONData = json::parse (stripComments (buffer.str ()));
        }
      catch (json::parse_error& e)
	{
	  cerr << "Unable to parse JSON file: " << e.what() << endl;
	  exit (1);
        }
      JSONFile.close ();
    }
  else
    {
      cerr << "Unable to open JSON file: " << filename << endl;
      exit (1);
    }

  return JSONData;
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

// ################################################################
// Function to check that directory exists, and create it otherwise
// ################################################################
bool Flux::CreateDirectory (const char* path)
{
  struct stat st = {0};
  
  if (stat (path, &st) == -1)
    {
#ifdef _WIN32
      if (mkdir (path) != 0)
	{
	  printf ("Error creating directory: %s\n", path);
	  return false;
	}
#else
      if (mkdir (path, 0700) != 0)
	{
	  printf ("Error creating directory: %s\n", path);
	  return false;
	}
#endif
    }
  
  return true;
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

