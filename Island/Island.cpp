#include "Island.h"

// ###########
// Constructor
// ###########
Island::Island ()
{
  // ----------------------------------------------
  // Ensure that directory ../Outputs/Island exists
  // ----------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Island"))
    {
      exit (1);
    }
  
  // ......................................
  // Read control parameters from JSON file
  // ......................................
  string JSONFilename = "../Inputs/Island.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  Nh   = JSONData["Nh"]  .get<int>    ();
  Nz   = JSONData["Nz"]  .get<int>    ();
  Nx   = JSONData["Nx"]  .get<int>    ();
  xmax = JSONData["xmax"].get<double> ();

  // ............
  // Sanity check
  // ............
  if (Nh < 1)
    {
      printf ("Island:: Error - Nh must be positive\n");
      exit (1);
    }
  if (Nz < 1)
    {
      printf ("Island:: Error - Nz must be positive\n");
      exit (1);
    }
  if (Nx < 1)
    {
      printf ("Island:: Error - Nx must be positive\n");
      exit (1);
    }
  if (xmax < 0.)
    {
      printf ("Island:: Error - xmax cannot be negative\n");
      exit (1);
    }
 
  printf ("\n");
  printf ("Class ISLAND::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Nh = %-4d Nz = %-4d Nx = %-4d xmax = %10.3e\n",
	  Nh, Nz, Nx, xmax);
}

// #########################
// Function to solve problem
// #########################
void Island::Solve ()
{
  // ..................
  // Set up radial grid
  // ..................
  xx = new double[Nx];
  for (int i = 0; i < Nx; i++)
    xx[i] = xmax * double (i) /double (Nx-1);

  // .............
  // Set up k grid
  // ..............
  kk = new double[Nx];
  double kmax = sqrt (4.*xmax*xmax + 1.) + 1.e-15;
  for (int i = 0; i < Nx; i++)
    kk[i] = 1. + (kmax - 1.) * double (i) /double (Nx-1);

  // .....................................................
  // Calculate perturbed temperature flux-surface function
  // .....................................................
  F = new double[Nx];

  double h, t_err;
  int    rept, count;

  double  k;
  double* Y   = new double[1];
  double* err = new double[1];

  rhs_chooser = 0;

  count = 0;
  h     = h0;
  k     = 1. + 1.e-15;
  Y[0]  = 0.;
  F[0]  = 0.;

  for (int i = 1; i < Nx; i++)
    {
      do
	{
	  CashKarp45Adaptive (1, k, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (k < kk[i]);
      CashKarp45Fixed (1, k, Y, err, kk[i] - k);

      F[i] = Y[0];
    }

  // .........................
  // Interpolate F(k) function
  // .........................
  Fspline = gsl_spline_alloc (gsl_interp_cspline, Nx);
  Facc    = gsl_interp_accel_alloc ();

  gsl_spline_init (Fspline, kk, F, Nx);

  // ...................
  // Calculate F_infinty
  // ...................
  double z;

  rhs_chooser = 1;
  
  count = 0;
  h     = h0;
  x     = xmax;
  Y[0]  = 0.;

  do
    {
      CashKarp45Adaptive (1, z, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (z < M_PI);
  CashKarp45Fixed (1, z, Y, err, M_PI - z);
  Finf = x - Y[0];

  printf ("F_infty = %11.4e\n", Finf);
   
  delete[] Y; delete[] err;

  // ............................................
  // Calculate harmonics of perturbed temperature
  // ............................................
  deltaTh.resize (Nh, Nx);

  double  zc;
  double* Y1   = new double[Nh];
  double* err1 = new double[Nh];

  rhs_chooser = 2;
  
  for (int j = 0; j < Nh; j++)
    deltaTh (j, 0) = 0.;

  for (int i = 1; i < Nx; i++)
    {
      count = 0;
      h     = h0;
      x     = xx[i];
      zc    = Getzetac (x);

      z = 0.;
      for (int j = 0; j < Nh; j++)
	Y1[j] = 0.;
 
      do
	{
	  CashKarp45Adaptive (Nh, z, Y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (z < zc);
      CashKarp45Fixed (Nh, z, Y1, err1, zc - z);

      for (int j = 0; j < Nh; j++)
	deltaTh (j, i) = Y1[j];
    }
  
  delete[] Y1; delete[] err1;

  // .......................
  // Verify small X behavior
  // .......................
  for (int j = 0; j < Nh; j++)
    {
      printf ("nu = %3d  (3/8) Y/X^3 = %11.4e\n", j, 3.*deltaTh (j, 1) /xx[1]/xx[1]/xx[1] /8.);
    }

  // ...................................................
  // Calculate temperature perturbation in z, zeta plane
  // ...................................................
  zz = new double[Nz];
  for (int j = 0; j < Nz; j++)
    zz[j] = 2.*M_PI * double (j) /double (Nz-1);

  dTo = new double[Nx];
  dTx = new double[Nx];
  deltaT.resize (Nx, Nz);
  
  for (int i = 0; i < Nx; i++)
    {
      double sum = 0., sum1 = 0.;

      for (int k = 0; k < Nh; k++)
	{
	  double nu = double (k);
	      
	  sum  += deltaTh (k, i) * cos (nu * M_PI);
	  sum1 += deltaTh (k, i);
	}
      dTo[i] = sum;
      dTx[i] = sum1;
      
      for (int j = 0; j < Nz; j++)
	{
	  double sum2 = 0.;
	  
	  for (int k = 0; k < Nh; k++)
	    {
	      double nu = double (k);
	      
	      sum2 += deltaTh (k, i) * cos (nu * zz[j]);
	    }
	  
	  deltaT (i, j) = sum2;
	}
    }
  
  // .........................
  // Write data to Netcdf file
  // .........................
  WriteNetcdf ();

  // ........
  // Clean up
  // ........
  delete[] xx; delete[] zz; delete[] dTo; delete[] kk; delete[] F; delete[] dTx;

  gsl_spline_free (Fspline);

  gsl_interp_accel_free (Facc);
}

// #####################################
// Function to write data to Netcdf file
// #####################################
void Island::WriteNetcdf ()
{
  printf ("Writing data to netcdf file Outputs/Layer/Layer.nc:\n");

   double* dTh_y = new double[Nh*Nx];
   double* dT_y  = new double[Nx*Nz];

   int cnt = 0;
   for (int i = 0; i < Nh; i++)
     for (int j = 0; j < Nx; j++)
       {
	 dTh_y[cnt] = deltaTh(i, j);
	 cnt++;
       }
   cnt = 0;
   for (int i = 0; i < Nx; i++)
     for (int j = 0; j < Nz; j++)
       {
	 dT_y[cnt] = deltaT(i, j);
	 cnt++;
       }
      
   try
     {
       NcFile dataFile ("../Outputs/Island/Island.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim h_d = dataFile.addDim ("Nh", Nh);
       NcDim x_d = dataFile.addDim ("Nx", Nx);
       NcDim z_d = dataFile.addDim ("Nz", Nz);
 
       vector<NcDim> T_d;
       T_d.push_back (h_d);
       T_d.push_back (x_d);

       vector<NcDim> TT_d;
       TT_d.push_back (x_d);
       TT_d.push_back (z_d);

       NcVar xx_x = dataFile.addVar  ("x",         ncDouble, x_d);
       xx_x.putVar (xx);
       NcVar kk_x = dataFile.addVar  ("k",         ncDouble, x_d);
       kk_x.putVar (kk);
       NcVar zz_x = dataFile.addVar  ("zeta",      ncDouble, z_d);
       zz_x.putVar (zz);
       NcVar F_x = dataFile.addVar   ("F",         ncDouble, x_d);
       F_x.putVar (F);
       NcVar dto_x = dataFile.addVar ("delta_T_o", ncDouble, z_d);
       dto_x.putVar (dTo);
       NcVar dtx_x = dataFile.addVar ("delta_T_x", ncDouble, z_d);
       dtx_x.putVar (dTx);
       NcVar dT_x = dataFile.addVar  ("delta_T_h", ncDouble, T_d);
       dT_x.putVar (dTh_y);
       NcVar dTT_x = dataFile.addVar ("delta_T",   ncDouble, TT_d);
       dTT_x.putVar (dT_y);
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/Island/Island.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }

   delete[] dTh_y; delete[] dT_y;
}

// #######################
// Function to calculate k
// #######################
double Island::Getk (double y, double zeta)
{
  double c = cos (zeta /2.);

  double kk = sqrt (c*c + 4.*y*y);

  if (kk >= 1. + 1.e-15)
    return kk;
  else
    return 1. + 1.e-15;
}

// ############################
// Function to calculate zeta_c
// ############################
double Island::Getzetac (double y)
{
  if (y < 0.5)
    return acos (1. - 8.*y*y);
  else
    return M_PI;
}

// ###################################
// Right hand sides of layer equations
// ###################################
void Island::CashKarp45Rhs (double z, double* Y, double* dYdz)
{
  if (rhs_chooser == 0)
    {
      double E = gsl_sf_ellint_Ecomp (1./z, GSL_MODE_DEFAULT);

      dYdz[0] = M_PI /4. /E;
    }
  else if (rhs_chooser == 1)
    {
      double k = Getk (x, z);
      double F = gsl_spline_eval (Fspline, k, Facc);

      dYdz[0] = F /M_PI;
    }
  else
    {
      double k = Getk (x, z);
      double E = gsl_sf_ellint_Ecomp (1./k, GSL_MODE_DEFAULT);
      double F = gsl_spline_eval (Fspline, k, Facc);

      dYdz[0] = F /M_PI;
      
      for (int nu = 1; nu < Nh; nu++)
	{
	  dYdz[nu] = (cos (double (nu-1) * z) - cos (double (nu+1) * z)) /16. /double (nu) /k /E;
	}
    }
}

