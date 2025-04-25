#include "Island.h"

#define NINPUT 6
#define NPARA 1

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

  Nh    = JSONData["Nh"]   .get<int>    ();
  Nb    = JSONData["Nb"]   .get<int>    ();
  Nz    = JSONData["Nz"]   .get<int>    ();
  NX    = JSONData["NX"]   .get<int>    ();
  Xmax  = JSONData["Xmax"] .get<double> ();
  delta = JSONData["delta"].get<double> ();

  // ............
  // Sanity check
  // ............
  if (Nh < 1)
    {
      printf ("Island:: Error - Nh must be positive\n");
      exit (1);
    }
  if (Nb < 1)
    {
      printf ("Island:: Error - Nb must be positive\n");
      exit (1);
    }
  if (Nz < 1)
    {
      printf ("Island:: Error - Nz must be positive\n");
      exit (1);
    }
  if (NX < 1)
    {
      printf ("Island:: Error - NX must be positive\n");
      exit (1);
    }
  if (Xmax < 1.)
    {
      printf ("Island:: Error - Xmax must be greater than unity\n");
      exit (1);
    }
  if (fabs(delta) > 1.)
    {
      printf ("Island:: Error - delta must lie in range -1 to +1\n");
      exit (1);
    }
 
  printf ("\n");
  printf ("Class ISLAND::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Nh = %-4d Nb = %-4d Nz = %-4d NX = %-4d Xmax = %10.3e delta = %10.3e\n",
	  Nh, Nb, Nz, NX, Xmax, delta);

  sprintf (buffer, "../Outputs/Island/Island.nc");
}

// #######################
// Alternative constructor
// #######################
Island::Island (int k, double _delta)
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

  Nh    = JSONData["Nh"]   .get<int>    ();
  Nb    = JSONData["Nb"]   .get<int>    ();
  Nz    = JSONData["Nz"]   .get<int>    ();
  NX    = JSONData["NX"]   .get<int>    ();
  Xmax  = JSONData["Xmax"] .get<double> ();

  delta = _delta;

  // ............
  // Sanity check
  // ............
  if (Nh < 1)
    {
      printf ("Island:: Error - Nh must be positive\n");
      exit (1);
    }
  if (Nb < 1)
    {
      printf ("Island:: Error - Nb must be positive\n");
      exit (1);
    }
  if (Nz < 1)
    {
      printf ("Island:: Error - Nz must be positive\n");
      exit (1);
    }
  if (NX < 1)
    {
      printf ("Island:: Error - NX must be positive\n");
      exit (1);
    }
  if (Xmax < 1.)
    {
      printf ("Island:: Error - Xmax must be greater than unity\n");
      exit (1);
    }
  if (fabs(delta) > 1.)
    {
      printf ("Island:: Error - delta must lie in range -1 to +1\n");
      exit (1);
    }
 
  printf ("\n");
  printf ("Class ISLAND::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("Rational surface %-4d: Nh = %-4d Nb = %-4d Nz = %-4d NX = %-4d Xmax = %10.3e delta = %10.3e\n",
	  k+1, Nh, Nb, Nz, NX, Xmax, delta);

  sprintf (buffer, "../Outputs/Island/Island%04d.nc", k);
}

// #########################
// Function to solve problem
// #########################
void Island::Solve (int FLAG)
{
  // .............
  // Set up k grid
  // ..............
  kk          = new double[NX];
  double Ymax = Xmax + delta /sqrt(8.);
  double kmax = sqrt (4.*Ymax*Ymax + 1.) + 1.e-15;
  for (int i = 0; i < NX; i++)
    kk[i] = 1. + (kmax - 1.) * double (i) /double (NX-1) + 1.e-15;

  // ..........................
  // Calculate E_n(k) functions
  // ..........................
  En.resize (Nb, NX);
  
  double h, t_err;
  int    rept, count;

  double  theta;
  double* Y   = new double[Nb];
  double* err = new double[Nb];

  rhs_chooser = 0;

  for (int i = 0; i < NX; i++)
    {
      X = 1./kk[i];

      count = 0;
      h     = h0;
      theta = 0.;

      for (int n = 0; n < Nb; n++)
	Y[n] = 0.;
      
      do
	{
	  CashKarp45Adaptive (Nb, theta, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (theta < M_PI/2.);
      CashKarp45Fixed (Nb, theta, Y, err, M_PI/2. - theta);

      for (int n = 0; n < Nb; n++)
	En(n, i) = Y[n];
    }
  delete[] Y; delete[] err;

  // ........................
  // Interpolate En functions
  // ........................
  Enspline = new gsl_spline* [Nb];
  Enacc    = new gsl_interp_accel* [Nb];

  for (int n = 0; n < Nb; n++)
    {
      Enspline[n] = gsl_spline_alloc (gsl_interp_cspline, NX);
      Enacc[n]    = gsl_interp_accel_alloc ();

      double* data = new double[NX];

      for (int i = 0; i < NX; i++)
	data[i] = En(n, i);

       gsl_spline_init (Enspline[n], kk, data, NX);

       delete[] data;
    }

  // ..............
  // Calculate dTdk
  // ..............
  dTdk = new double[NX];

  for (int i = 0; i < NX; i++)
    dTdk[i] = (M_PI/4.) /GetG (kk[i]);

  // ................
  // Interpolate dTdk
  // ................
  dTdkspline = gsl_spline_alloc (gsl_interp_cspline, NX);
  dTdkacc    = gsl_interp_accel_alloc ();

  gsl_spline_init (dTdkspline, kk, dTdk, NX);

  // ..................
  // Set up radial grid
  // ..................
  XX = new double[NX];
  for (int i = 0; i < NX; i++)
    XX[i] = - Xmax + 2.*Xmax * double (i) /double (NX-1);

  // .....................................................
  // Calculate perturbed temperature flux-surface function
  // .....................................................
  F = new double[NX];

  double  k;
  Y   = new double[1];
  err = new double[1];

  rhs_chooser = 1;

  count = 0;
  h     = h0;
  k     = 1. + 1.e-15;
  Y[0]  = 0.;
  F[0]  = 0.;

  for (int i = 1; i < NX; i++)
    {
      do
	{
	  CashKarp45Adaptive (1, k, Y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (k < kk[i]);
      CashKarp45Fixed (1, k, Y, err, kk[i] - k);

      F[i] = Y[0];
    }
   delete[] Y; delete[] err;

  // .........................
  // Interpolate F(k) function
  // .........................
  Fspline = gsl_spline_alloc (gsl_interp_cspline, NX);
  Facc    = gsl_interp_accel_alloc ();

  gsl_spline_init (Fspline, kk, F, NX);

  // ............................................
  // Calculate harmonics of perturbed temperature
  // ............................................
  deltaTh.resize (Nh, NX);

  double  xi, xic;
  double* Y1   = new double[Nh];
  double* err1 = new double[Nh];
          ximx = new double[NX];

  rhs_chooser = 2;
  
  for (int j = 0; j < Nh; j++)
    deltaTh (j, 0) = 0.;

  for (int i = 0; i < NX; i++)
    {
      if (i%100 == 0)
	{
	  printf (".");
	  fflush (stdout);
	}
      
      count   = 0;
      h       = h0;
      X       = XX[i];
      xic     = GetXic (X);
      ximx[i] = xic /M_PI;

      xi = 0.;
      for (int j = 0; j < Nh; j++)
	Y1[j] = 0.;
 
      do
	{
	  CashKarp45Adaptive (Nh, xi, Y1, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	}
      while (xi < xic);
      CashKarp45Fixed (Nh, xi, Y1, err1, xic - xi);

      for (int j = 0; j < Nh; j++)
	deltaTh (j, i) = Y1[j];
    }
  delete[] Y1; delete[] err1;

  // ....................
  // Calculate T0_infinty
  // ....................
  T0inf = XX[NX-1] - deltaTh (0, NX-1);
  
  printf (" T0_infty = %11.4e\n", T0inf);

  if (FLAG)
    {
      // ...................................................
      // Calculate temperature perturbation in z, zeta plane
      // ...................................................
      zz = new double[Nz];
      for (int j = 0; j < Nz; j++)
	zz[j] = 2.*M_PI * double (j) /double (Nz-1);
      
      dTo = new double[NX];
      dTx = new double[NX];
      deltaT.resize (NX, Nz);
      
      for (int i = 0; i < NX; i++)
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
    }
      
  // .........................
  // Write data to Netcdf file
  // .........................
  WriteNetcdf (FLAG);

  // ........
  // Clean up
  // ........
  printf ("Cleaning up\n");
  
  delete[] XX; delete[] kk; delete[] dTdk; delete[] ximx; delete[] F; 

  if (FLAG)
    {
       delete[] zz;
       delete[] dTo;
       delete[] dTx;
    }

  for (int n = 0; n < Nb; n++)
    {
      gsl_spline_free (Enspline[n]);

      gsl_interp_accel_free (Enacc[n]);
    }
  delete[] Enspline;
  delete[] Enacc;

  gsl_spline_free (dTdkspline);
  gsl_spline_free (Fspline);

  gsl_interp_accel_free (dTdkacc);
  gsl_interp_accel_free (Facc);
}

// #####################################
// Function to write data to Netcdf file
// #####################################
void Island::WriteNetcdf (int FLAG)
{
  printf ("Writing data to netcdf file %s\n", buffer);

  double Input[NINPUT], Para[NPARA];

  double* En_y  = new double[Nb*NX];
  double* dTh_y = new double[Nh*NX];
  double* dT_y  = new double[NX*Nz];

  Input[0] = double (Nh);
  Input[1] = double (Nb);
  Input[2] = double (Nz);
  Input[3] = double (NX);
  Input[4] = double (Xmax);
  Input[5] = double (delta);
   
  Para[0] = T0inf;

  int cnt = 0;
  for (int i = 0; i < Nb; i++)
    for (int j = 0; j < NX; j++)
      {
	En_y[cnt] = En(i, j);
	cnt++;
      }
  
  cnt = 0;
  for (int i = 0; i < Nh; i++)
    for (int j = 0; j < NX; j++)
      {
	dTh_y[cnt] = deltaTh(i, j);
	cnt++;
      }

  if (FLAG)
    {
      cnt = 0;
      for (int i = 0; i < NX; i++)
	for (int j = 0; j < Nz; j++)
	  {
	    dT_y[cnt] = deltaT(i, j);
	    cnt++;
	  }
    }

   try
    {
      NcFile dataFile (buffer, NcFile::replace);
 
      dataFile.putAtt ("Git_Hash",     GIT_HASH);
      dataFile.putAtt ("Compile_Time", COMPILE_TIME);
      dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

      NcDim i_d = dataFile.addDim ("Ni", NINPUT);
      NcDim p_d = dataFile.addDim ("Np", NPARA);
      NcDim h_d = dataFile.addDim ("Nh", Nh);
      NcDim b_d = dataFile.addDim ("Nb", Nh);
      NcDim x_d = dataFile.addDim ("NX", NX);
      NcDim z_d = dataFile.addDim ("Nz", Nz);

      vector<NcDim> E_d;
      E_d.push_back (b_d);
      E_d.push_back (x_d);
            
      vector<NcDim> T_d;
      T_d.push_back (h_d);
      T_d.push_back (x_d);
      
      vector<NcDim> TT_d;
      TT_d.push_back (x_d);
      TT_d.push_back (z_d);

      NcVar i_x    = dataFile.addVar ("InputParameters", ncDouble, i_d);
      i_x.putVar (Input);
      NcVar p_x    = dataFile.addVar ("para",            ncDouble, p_d);
      p_x.putVar (Para);
      NcVar xx_x   = dataFile.addVar ("X",               ncDouble, x_d);
      xx_x.putVar (XX);
      NcVar kk_x   = dataFile.addVar ("k",               ncDouble, x_d);
      kk_x.putVar (kk);
      NcVar En_x   = dataFile.addVar ("E_n",             ncDouble, E_d);
      En_x.putVar (En_y);
      NcVar dTdk_x = dataFile.addVar ("dTdk",            ncDouble, x_d);
      dTdk_x.putVar (dTdk);
      NcVar F_x    = dataFile.addVar ("F",               ncDouble, x_d);
      F_x.putVar (F);
      NcVar dT_x   = dataFile.addVar ("delta_T_h",       ncDouble, T_d);
      dT_x.putVar (dTh_y);
      NcVar xim_x  = dataFile.addVar ("xi_c",            ncDouble, x_d);
      xim_x.putVar (ximx);

      if (FLAG)
	{
	  NcVar zz_x  = dataFile.addVar ("zeta",      ncDouble, z_d);
	  zz_x.putVar (zz);
 	  NcVar dto_x = dataFile.addVar ("delta_T_o", ncDouble, x_d);
	  dto_x.putVar (dTo);
	  NcVar dtx_x = dataFile.addVar ("delta_T_x", ncDouble, x_d);
	  dtx_x.putVar (dTx);	
	  NcVar dTT_x = dataFile.addVar ("delta_T",   ncDouble, TT_d);
	  dTT_x.putVar (dT_y);
	}
 
      dataFile.close ();
     }
  catch (NcException& e)
    {
      printf ("Error writing data to netcdf file %s\n", buffer);
      printf ("%s\n", e.what ());
      exit (1);
    }

  delete[] En_y; delete[] dTh_y; delete[] dT_y;
}

// ##########################
// Function to calculate G(k)
// ##########################
double Island::GetG (double k)
{
  double G = gsl_spline_eval (Enspline[0], k, Enacc[0]);

  for (int n = 1; n < Nb; n++)
    {
      double nn = double (n);

      G += 2. * cos (nn*M_PI) * gsl_sf_bessel_Jn (n, nn*delta*delta) * gsl_spline_eval (Enspline[n], k, Enacc[n]);
    }

  return G;
}

// #############################
// Function to return cos (zeta)
// #############################
double Island::GetCosZeta (double xi)
{
  double sum = - delta*delta/2.;

  for (int j = 1; j < Nb; j++)
    {
      double jj = double (j);

      sum += cos (jj * xi) * (gsl_sf_bessel_Jn (j-1, jj*delta*delta) - gsl_sf_bessel_Jn (j+1, jj*delta*delta)) /jj;
    }

  return sum;
}

// #############################
// Function to return sin (zeta)
// #############################
double Island::GetSinZeta (double xi)
{
  if (fabs (delta) < 1.e-15)
    return sin (xi);
  else
    {
      double sum = 0.;

      for (int j = 1; j < Nb; j++)
	{
	  double jj = double (j);
	  
	  sum += 2. * sin (jj * xi) * gsl_sf_bessel_Jn (j, jj*delta*delta) /jj /delta/delta;
	}
   
      return sum;
    }
}

// ###############################
// Function to return cos (n*zeta)
// ###############################
double Island::GetCosnZeta (int n, double xi)
{
  double nn  = double (n);
  double sum = 0.;

  for (int j = 1; j < Nb; j++)
    {
      double jj = double (j);

      if (j - n >= 0)
	sum += nn * cos (jj * xi) * (gsl_sf_bessel_Jn (j-n, jj*delta*delta) - gsl_sf_bessel_Jn (j+n, jj*delta*delta)) /jj;
      else
	sum += nn * cos (jj * xi) * (cos (M_PI * double (j-n)) * gsl_sf_bessel_Jn (n-j, jj*delta*delta) - gsl_sf_bessel_Jn (j+n, jj*delta*delta)) /jj;
    }

  return sum;
}

// ##############################
// Function to return sin(n*zeta)
// ##############################
double Island::GetSinnZeta (int n, double xi)
{
  double nn  = double (n);
  double sum = 0.;

  for (int j = 1; j < Nb; j++)
    {
      double jj = double (j);

      if (j - n >= 0)
	sum += nn * sin (jj * xi) * (gsl_sf_bessel_Jn (j-n, jj*delta*delta) + gsl_sf_bessel_Jn (j+n, jj*delta*delta)) /jj;
      else
	sum += nn * sin (jj * xi) * (cos (M_PI * double (j-n)) * gsl_sf_bessel_Jn (n-j, jj*delta*delta) + gsl_sf_bessel_Jn (j+n, jj*delta*delta)) /jj;
    }

  return sum;
}

// ###############################
// Function to calculate xi (zeta)
// ###############################
double Island::GetXi (double zeta)
{
  return zeta - delta*delta * sin (zeta);
}

// ###############################
// Function to calculate zeta (xi)
// ###############################
double Island::GetZeta (double xi)
{
  double sum = xi;

  for (int j = 1; j < Nb; j++)
    {
      double jj = double (j);
      
      sum += 2. * sin (jj * xi) * gsl_sf_bessel_Jn (j, jj*delta*delta) /jj;
    }

  return sum;
}

// ##############################
// Function to calculate xi_c (X)
// ##############################
double Island::GetXic (double x)
{
  if (x < - 0.5 - delta/sqrt(8.) || x > 0.5 - delta/sqrt(8.))
    return M_PI;
  else
    {
      X = x;
      
      double z = RootFind (0., M_PI);

      return z;
    }
}
  
// ###############################
// Function to calculate k (X, xi)
// ###############################
double Island::Getk (double x, double xi)
{
  double zeta = GetZeta (xi);
  double Y    = x - (delta /sqrt(8.)) * cos (zeta);
  double c    = cos (xi /2.);
  double kx   = sqrt (c*c + 4.*Y*Y);

  if (kx < kk[0])
    return kk[0];
  else if (kx > kk[NX-1])
    return kk[NX-1];
  else
    return kx;
}

// ################################
// Function to calculate sigma (xi)
// ################################
double Island::GetSigma (double xi)
{
  double sum = 1.;

  for (int j = 1; j < Nb; j++)
    {
      double jj = double (j);

      sum += 2. * cos (jj * xi) * gsl_sf_bessel_Jn (j, jj*delta*delta);
    }

  return sum;
}

// ###################################
// Function to calculate kappa (X, xi)
// ###################################
double Island::GetKappa (double x, double xi)
{
  double zeta = GetZeta (xi);
  
  return sin (xi) * (1. - delta*delta * cos (zeta)) - 2.*sqrt(8.) * delta * x * sin (zeta) + delta*delta * sin (2.*zeta);
}

// ###################################
// Right hand sides of layer equations
// ###################################
void Island::CashKarp45Rhs (double z, double* Y, double* dYdz)
{
  if (rhs_chooser == 0)
    {
      for (int n = 0; n < Nb; n++)
	{
	  double nn = double (n);

	  dYdz[n] = cos (2.*nn*z) * sqrt (1. - X*X * sin (z) * sin (z));
	}
    }
  else if (rhs_chooser == 1)
    {
      if (z > kk[NX-1])
	dYdz[0] = dTdk[NX-1];
      else
	dYdz[0] = gsl_spline_eval (dTdkspline, z, dTdkacc);
    }
  else
    {
      double k     = Getk (X, z);
      double sigma = GetSigma (z);
      double kappa = GetKappa (X, z);
      double G     = GetG (k);
      double F     = gsl_spline_eval (Fspline, k, Facc);
      double Y     = X - delta * GetCosZeta (z) /sqrt(8.);

      if (Y > 0.)
	dYdz[0] =   sigma * F /M_PI;
      else
	dYdz[0] = - sigma * F /M_PI;
      
      for (int nu = 1; nu < Nh; nu++)
	{
	  double sinnu;
	  if (nu == 1)
	    sinnu = GetSinZeta (z);
	  else
	    sinnu = GetSinnZeta (nu, z);

	  if (Y > 0.)
	    dYdz[nu] =   sinnu * kappa * sigma /k /G /8. /double (nu);
	  else
	    dYdz[nu] = - sinnu * kappa * sigma /k /G /8. /double (nu);
	}
    }
}

// ################################
// Target function for finding xi_c
// ################################
double Island::RootFindF (double xi)
{ 
  double cosz = GetCosZeta (xi);

  return 1. - (sqrt(8.)*X - delta * cosz) * (sqrt(8.)*X - delta * cosz) - cos (xi);
}
