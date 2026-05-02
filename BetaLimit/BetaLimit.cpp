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

  CNTRL     = JSONData["CNTRL"]    .get<int>();
  UP        = JSONData["UP"]       .get<int>();
  val_start = JSONData["val_start"].get<double>();
  val_end   = JSONData["val_end"]  .get<double>();
  N_val     = JSONData["N_val"]    .get<int>();
  
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
  double pc_en = pc_end;
  for (int i = 0; i <= N_val; i++)
    { 
      val = val_start + (val_end - val_start) * double (i) /double (N_val);
 
      double pc = RootFind (pc_st, pc_en);

      if (UP)
	{
	  pc_st = pc;
	  pc_en = pc_end;
	}
      else
	{
	  pc_st = pc;
	  pc_en = pc_start;
	}

      fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", val, rs, ss, pc, dW);
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

  tj.Solve (CNTRL, val, _pc);

  try
    {
      NcFile dataFile ("../Outputs/TJ/TJ.nc", NcFile::read);

      NcVar W_x   = dataFile.getVar ("delta_W");
      NcVar Wv_x  = dataFile.getVar ("delta_W_v");
      NcVar pW_x  = dataFile.getVar ("pdelta_W");
      NcVar pWv_x = dataFile.getVar ("pdelta_W_v");
      NcDim W_d   = W_x.getDim (0);

      int     NW       = W_d.getSize ();
      double* deltaW   = new double[NW];
      double* deltaWv  = new double[NW];
      double* pdeltaW  = new double[NW];
      double* pdeltaWv = new double[NW];

      W_x  .getVar (deltaW);
      Wv_x .getVar (deltaWv);
      pW_x .getVar (pdeltaW);
      pWv_x.getVar (pdeltaWv);

      if (CNTRL == 1)
	{
	  int i = -1;
	  do
	    {
	      i++;
	      dW = pdeltaW[i];
	    }
	  while (pdeltaWv[i] < 0.);
	}
      else
	{
	  int i = -1;
	  do
	    {
	      i++;
	      dW = deltaW[i];
	    }
	  while (deltaWv[i] < 0.);
	}
  
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

      delete[] deltaW; delete[] deltaWv; delete[] pdeltaW; delete[] pdeltaWv;
      delete[] r_res; delete[] s_res;
    }
  catch (NcException& e)
    {
      printf ("Error reading data from netcdf file Outputs/TJ/TJ.nc\n");
      printf ("%s\n", e.what ());
      exit (1);
    }

  FILE* file = OpenFilea ("../Outputs/BetaLimit/Record.out");

  fprintf (file, "%11.4e %11.4e %11.4e %11.4e %11.4e\n", val, _pc, rs, ss, dW);

  fclose (file);

  return dW;
}
