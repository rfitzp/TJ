// BetaLimit.cpp

#include "BetaLimit.h"

// ###########
// Constructor
// ###########
BetaLimit::BetaLimit ()
{
  // --------------------------------------
  // Read control parameters from JSON file
  // --------------------------------------
  string JSONFilename = "../Inputs/BetaLimit.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  qc_start = JSONData["qc_start"].get<double>();
  qc_end   = JSONData["qc_end"]  .get<double>();
  N_qc     = JSONData["N_qc"]    .get<int>();
  
  pc_start = JSONData["pc_start"].get<double>();
  pc_end   = JSONData["pc_end"]  .get<double>();
  Nint     = JSONData["Nint"]    .get<int>();

  Eta     = 1.e-6;
  Maxiter = 20;
}

// ##########
// Destructor
// ##########
BetaLimit::~BetaLimit ()
{
}

// #########################
// Function to solve problem
// #########################
void BetaLimit::Solve ()
{
  FILE* file = OpenFilew ("../Outputs/BetaLimit/BetaLimit.out");

  double pc_st = pc_start;
  for (int i = 0; i < N_qc; i++)
    { 
      qc = qc_start + (qc_end - qc_start) * double (i) /double (N_qc - 1);

      double pc = RootFind (pc_st, pc_end);
      
      pc_st = pc;

      fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", qc, rs, ss, pc, dW);
      fflush (file);
    }
  
  fclose (file);
}

// ################################
// Target function for root finding
// ################################
double BetaLimit::RootFindF (double _pc)
{
  TJ tj;

  tj.Solve (qc, _pc);

  try
    {
      NcFile dataFile ("../Outputs/TJ/TJ.nc", NcFile::read);

      NcVar W_x = dataFile.getVar ("delta_W");
      NcDim W_d = W_x.getDim (0);

      int     NW     = W_d.getSize ();
      double* deltaW = new double[NW];

      W_x.getVar (deltaW);

      dW = deltaW[1];

      NcVar r_x = dataFile.getVar ("r_res");
      NcDim r_d = r_x.getDim (0);

      int     Nr     = r_d.getSize ();
      double* r_res = new double[Nr];

      r_x.getVar (r_res);

      rs = r_res[0];

      NcVar s_x = dataFile.getVar ("s_res");
      NcDim s_d = s_x.getDim (0);

      int     Ns    = s_d.getSize ();
      double* s_res = new double[Ns];

      s_x.getVar (s_res);

      ss = s_res[0];

      delete [] deltaW; delete[] r_res; delete[] s_res;
    }
  catch (NcException& e)
    {
      printf ("Error reading data from netcdf file Outputs/TJ/TJ.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  FILE* file = OpenFilea ("../Outputs/BetaLimit/Record.out");

  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", qc, _pc, rs, ss, dW);

  fclose (file);

  return dW;
}
