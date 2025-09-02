#include "FourField.h"

#define NINPUT 15

// ###########
// Constructor
// ###########
FourField::FourField ()
{
  // -------------------------------------------------
  // Ensure that directory ../Outputs/FourField exists
  // -------------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/FourField"))
    {
      exit (1);
    }

  // ......................................
  // Read control parameters from JSON file
  // ......................................
  string JSONFilename = "../Inputs/FourField.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  g_r   = JSONData["g_r"]  .get<double> ();
  g_i   = JSONData["g_i"]  .get<double> ();
  Qe    = JSONData["Qe"]   .get<double> ();
  Qi    = JSONData["Qi"]   .get<double> ();
  D     = JSONData["D"]    .get<double> ();
  Pphi  = JSONData["Pphi"] .get<double> ();
  Pperp = JSONData["Pperp"].get<double> ();
  cbeta = JSONData["cbeta"].get<double> ();

  JK = JSONData["JK"].get<int> ();

  if (JK)
    Pperp = cbeta*cbeta;

  iotae = - Qe /(Qi - Qe);

  pstart = JSONData["pstart"].get<double> ();
  pend   = JSONData["pend"]  .get<double> ();
  P3max  = JSONData["P3max"] .get<double> ();

  Scan  = JSONData["Scan"] .get<int>    ();
  Fscan = JSONData["Fscan"].get<double> ();
  Nscan = JSONData["Nscan"].get<int>    ();

  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  // ------------
  // Sanity check
  // ------------
  if (D < 0.)
    {
      printf ("FourField:: Error - D cannot be negative\n");
      exit (1);
    }
  if (Pphi < 0.)
    {
      printf ("FourField:: Error - Pphi cannot be negative\n");
      exit (1);
    }
  if (Pperp < 0.)
    {
      printf ("FourField:: Error - Pperp cannot be negative\n");
      exit (1);
    }
  if (cbeta < 0.)
    {
      printf ("FourField:: Error - cbeta cannot be negative\n");
      exit (1);
    }
  if (pend < 0.)
    {
      printf ("FourField:: Error - pstart cannot be negative\n");
      exit (1);
    }
  if (pstart < pend)
    {
      printf ("FourField:: Error - pstart cannot be less than pend\n");
      exit (1);
    }
  if (P3max < 0.)
    {
      printf ("FourField:: Error - P3max cannot be negative\n");
      exit (1);
    }
  if (acc < 0.)
    {
      printf ("FourField:: Error - acc cannot be negative\n");
      exit (1);
    }
  if (h0 < 0.)
    {
      printf ("FourField:: Error - h0 cannot be negative\n");
      exit (1);
    }
   if (hmin < 0.)
    {
      printf ("FourField:: Error - hmin cannot be negative\n");
      exit (1);
    }
   if (hmax < hmin)
    {
      printf ("FourField:: Error - hmax cannot be less than hmin\n");
      exit (1);
    }

  printf ("\n");
  printf ("Class FOURFIELD::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("g_r    = %10.3e g_i   = %10.3e Qe    = %10.3e Qi    = %10.3e D = %10.3e JK = %1d\n",
	  g_r, g_i, Qe, Qi, D, JK);
  printf( "Pphi   = %10.3e Pperp = %10.3e cbeta = %10.3e iotae = %10.3e\n",
	  Pphi, Pperp, cbeta, iotae);
  printf ("pstart = %10.3e pend  = %10.3e P3max = %10.3e\n",
	  pstart, pend, P3max);
  printf ("acc    = %10.3e h0    = %10.3e hmin  = %10.3e hmax  = %10.3e\n",
	  acc, h0, hmin, hmax);
}

// #########################
// Function to solve problem
// #########################
void FourField::Solve ()
{
  // .................................
  // Solve three-field layer equations
  // .................................
  SolveThreeFieldLayerEquations ();
  printf ("\nDeltas3 = (%10.3e, %10.3e)\n", real(Deltas3), imag(Deltas3));

  // ................................
  // Solve four-field layer equations
  // ................................
  SolveFourFieldLayerEquations ();
  printf ("Deltas4 = (%10.3e, %10.3e)\n", real(Deltas4), imag(Deltas4));

  // ......................
  // Perform frequency scan
  // ......................
  if (Scan)
    {
      printf ("\n");
      
      givals  = new double[Nscan+1];
      D3rvals = new double[Nscan+1];
      D3ivals = new double[Nscan+1];
      D4rvals = new double[Nscan+1];
      D4ivals = new double[Nscan+1];
      J1rvals = new double[Nscan+1];
      J1ivals = new double[Nscan+1];
      J2rvals = new double[Nscan+1];
      J2ivals = new double[Nscan+1];

      double gi_save = g_i;
      double gi_beg  = - Fscan * (Qe - Qi);
      double gi_end  = + Fscan * (Qe - Qi);

      complex<double> _DeltaJ1, _DeltaJ2;

      for (int i = 0; i <= Nscan; i++)
	{
	  g_i = gi_beg + double (i) * (gi_end - gi_beg) /double (Nscan);

	  SolveThreeFieldLayerEquations ();
	  SolveFourFieldLayerEquations  ();

	  GetJace (_DeltaJ1, _DeltaJ2);

	  givals [i] = g_i;
	  D3rvals[i] = real(Deltas3);
	  D3ivals[i] = imag(Deltas3);
	  D4rvals[i] = real(Deltas4);
	  D4ivals[i] = imag(Deltas4);
	  J1rvals[i] = real(_DeltaJ1);
	  J1ivals[i] = imag(_DeltaJ1);
	  J2rvals[i] = real(_DeltaJ2);
	  J2ivals[i] = imag(_DeltaJ2);

	  if (i % 10 == 0)
	    printf ("%4d g = (%10.3e, %10.3e) Delta3 = (%10.3e, %10.3e) Delta4 = (%10.3e, %10.3e)\n",
		  i, g_r, givals[i], D3rvals[i], D3ivals[i], D4rvals[i], D4ivals[i]);
	}

      g_i = gi_save;
    }
     
  // .....................................
  // Write calculation date to netcdf file
  // .....................................
  WriteNetcdf ();
  
  // ........
  // Clean up
  // ........
  CleanUp ();
}

// ######################################
// Function to write  data to netcdf file
// ######################################
void FourField::WriteNetcdf ()
{
   printf ("Writing data to netcdf file Outputs/FourField/FourField.nc:\n");

   double Input[NINPUT];
   
   Input[0]  = g_r;
   Input[1]  = g_i;
   Input[2]  = Qe;
   Input[3]  = Qi;
   Input[4]  = D;
   Input[5]  = Pphi;
   Input[6]  = Pperp;
   Input[7]  = cbeta;
   Input[8]  = pstart;
   Input[9]  = pend;
   Input[10] = P3max;
   Input[11] = acc;
   Input[12] = h0;
   Input[13] = hmin;
   Input[14] = hmax;

   try
     {
       NcFile dataFile ("../Outputs/FourField/FourField.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim i_d = dataFile.addDim ("Ni",    NINPUT);
       NcDim x_d = dataFile.addDim ("Nscan", Nscan+1);

       NcVar i_x = dataFile.addVar ("InputParameters", ncDouble, i_d);
       i_x.putVar (Input);

       if (Scan)
	 {
	   NcVar givals_x  = dataFile.addVar ("g_i",             ncDouble, x_d);
	   givals_x.putVar (givals);
	   NcVar D3rvals_x = dataFile.addVar ("Delta3_r",        ncDouble, x_d);
	   D3rvals_x.putVar (D3rvals);
	   NcVar D3ivals_x = dataFile.addVar ("Delta3_i",        ncDouble, x_d);
	   D3ivals_x.putVar (D3ivals);
	   NcVar D4rvals_x = dataFile.addVar ("Delta4_r",        ncDouble, x_d);
	   D4rvals_x.putVar (D4rvals);
	   NcVar D4ivals_x = dataFile.addVar ("Delta4_i",        ncDouble, x_d);
	   D4ivals_x.putVar (D4ivals);
	   NcVar J1rvals_x = dataFile.addVar ("DeltaJ1_r",       ncDouble, x_d);
	   J1rvals_x.putVar (J1rvals);
	   NcVar J1ivals_x = dataFile.addVar ("DeltaJ1_i",       ncDouble, x_d);
	   J1ivals_x.putVar (J1ivals);
	   NcVar J2rvals_x = dataFile.addVar ("DeltaJ2_r",       ncDouble, x_d);
	   J2rvals_x.putVar (J2rvals);
	   NcVar J2ivals_x = dataFile.addVar ("DeltaJ2_i",       ncDouble, x_d);
	   J2ivals_x.putVar (J2ivals);
	 }
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/FourField/FourField.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }
}

// #############################
// Function to deallocate memory
// #############################
void FourField::CleanUp ()
{
  printf ("Cleaning up\n");

  if (Scan)
    {
      delete[] givals;

      delete[] D3rvals; delete[] D3ivals; delete[] D4rvals; delete[] D4ivals;
      delete[] J1rvals; delete[] J1ivals; delete[] J2rvals; delete[] J2ivals;
    }
}

// #############################################
// Function to solve three-field layer equations
// #############################################
void FourField::SolveThreeFieldLayerEquations ()
{
  // ........
  // Define g
  // ........
  complex<double> g = complex<double> (g_r, g_i);

  // ..............
  // Determine Pmax
  // ..............
  complex<double> gEe = g + mc_i * Qe;
  complex<double> gEi = g + mc_i * Qi;
  complex<double> gPD = Pperp + (g + mc_i * Qi) * D*D;
  double          PS  = Pphi + Pperp;
  double          PP  = Pphi * Pperp;
  double          PD  = Pphi * D*D /iotae;
  
  double Pmax[6];
  Pmax[0] = pow (abs (gEe),          0.5);
  Pmax[1] = pow (abs (gEi * PS /PP), 0.5);
  Pmax[2] = pow (abs (g * gEi /PP),  0.25);
  Pmax[3] = pow (abs (gPD /PD),      0.5);
  Pmax[4] = pow (abs (gEe /PD),      0.25);
  Pmax[5] = pow (abs (PD /PP),       0.25);

  if (Pmax[3] < P3max)
    lowD = 0;
  else
    {
      lowD    = 1;
      Pmax[3] = pow (abs (gEe /Pperp), 0.5);
      Pmax[4] = pow (Pphi, -1./6.);
      Pmax[5] = 1.;
    }
 
  double PMAX = 1.;
  for (int i = 0; i < 6; i++)
    if (Pmax[i] > PMAX)
      PMAX = Pmax[i];

  // ...........................................................
  // Integrate layer equations backward in p to calculate Deltas
  // ...........................................................
  complex<double>   alpha, beta, gamma, X;
  double            p, h, t_err;
  int               rept; count = 0; rhs_chooser = 0;
  complex<double>*  y    = new complex<double>[1];
  complex<double>*  dydp = new complex<double>[1];
  complex<double>*  err  = new complex<double>[1];

  alpha = - gEe;
  p     = pstart * PMAX;
  h     = - h0;
  if (!lowD)
    {
      beta  = PP /PD;
      gamma = beta * (1. + gEi * PS /PP - gPD /PD);
      X     = (gamma - sqrt (beta) * (1. - sqrt (beta) * alpha)) / (2. * sqrt (beta));
      
      y[0] = X - sqrt (beta) * p*p;
    }
  else
    {
      beta  = Pphi;
      gamma = - mc_i * (Qe - Qi) * Pphi /Pperp + gEi;
      X     = (alpha * beta - gamma) /2./ sqrt(beta);

      y[0] = - 1. + X * p - sqrt(beta) * p*p*p;
    }
  
  do
    {
      CashKarp45Adaptive (1, p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (p > pend);
  CashKarp45Fixed (1, p, y, err, pend - p);
  
  CashKarp45Rhs (p, y, dydp);
  Deltas3 = M_PI /dydp[0];

  // ........
  // Clean up
  // ........
  delete[] y; delete[] dydp; 
}

// #############################################
// Function to solve four-field layer equations
// #############################################
void FourField::SolveFourFieldLayerEquations ()
{
  // ........
  // Define g
  // ........
  complex<double> g = complex<double> (g_r, g_i);

  // ..............
  // Determine Pmax
  // ..............
  complex<double> gEe  = g + mc_i * Qe;
  complex<double> gEi  = g + mc_i * Qi;
  complex<double> gPD  = Pperp + (g + mc_i * Qi) * D*D;
  complex<double> R    = Pperp + (1. - 1./iotae) * D*D * gEi;
  double          PD   = Pphi * D*D /iotae;
  double          cbm2 = 1./cbeta/cbeta;

  complex<double> F11_4 = gEi + Pphi * gEe;
  complex<double> F12_4 = gEi + Pphi * gEe /iotae;
  complex<double> F21_4 = - mc_i * Qe * Pphi        + cbm2 * (g * D*D * gEi + Pphi * mc_i * Qe);
  complex<double> F22_4 = - mc_i * Qe * Pphi /iotae + cbm2 * (g * gPD       + Pphi * gEe);

  complex<double> F11_6 = Pphi;
  complex<double> F12_6 = Pphi /iotae;
  complex<double> F21_6 = cbm2 * D*D * g * Pphi        + cbm2 * D*D * gEi * Pphi;
  complex<double> F22_6 = cbm2 * D*D * g * Pphi /iotae + cbm2 * gPD       * Pphi;

  complex<double> F21_8 = cbm2 * D*D * Pphi*Pphi;
  complex<double> F22_8 = cbm2 * D*D * Pphi*Pphi /iotae;

  double Pmax[4];
  Pmax[0] = pow (abs (F21_6/F21_8), 0.5);
  Pmax[1] = pow (abs (F22_6/F22_8), 0.5);
  Pmax[2] = pow (abs (F11_4/F11_6), 0.5);
  Pmax[3] = pow (abs (F12_4/F12_6), 0.5);
  //Pmax[4] = pow (abs (F21_4/F21_6), 0.5);
  //Pmax[5] = pow (abs (F22_4/F22_6), 0.5);

  double PMAX = 1.;
  for (int i = 0; i < 4; i++)
    if (Pmax[i] > PMAX)
      PMAX = Pmax[i];

  // ...........................................................
  // Integrate layer equations backward in p to calculate Deltas
  // ...........................................................
  double            p, h, t_err;
  int               rept; count = 0; rhs_chooser = 1;
  complex<double>*  y    = new complex<double>[4];
  complex<double>*  dydp = new complex<double>[4];
  complex<double>*  err  = new complex<double>[4];
 
  p     = pstart * PMAX;
  h     = - h0;

  y[0] = - sqrt (iotae) * sqrt (R) * p*p /D - cbeta * sqrt (iotae) * p*p /D;
  y[1] = - cbeta * p*p /D /sqrt (iotae);
  y[2] = - sqrt (iotae) * D * Pphi * p*p*p*p /cbeta;
  y[3] = - D * Pphi * p*p*p*p /sqrt (iotae) /cbeta;
  
  do
    {
      CashKarp45Adaptive (4, p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (p > pend);
  CashKarp45Fixed (4, p, y, err, pend - p);
 
  CashKarp45Rhs (p, y, dydp);
  Deltas4 = M_PI /(dydp[0] - dydp[1] * mc_i * Qe /gEe);

  /*
  complex<double> E12 = - 2. * mc_i * Qe /(complex<double> (g_r, g_i) + mc_i * Qe);
  printf ("W_11 = (%10.3e, %10.3e) W_21 = (%10.3e, %10.3e)\n",
	  real (y[0]), imag(y[0]), real (2.*y[2]/E12), imag (2.*y[2]/E12));
  complex<double> W21 = 0.5 * E12 * (dydp[0] - dydp[3]);
  printf ("W_21 = (%10.3e, %10.3e) W_21 = (%10.3e, %10.3e)\n",
	  real (dydp[2]), imag(dydp[2]), real (W21), imag (W21));
  */

  // ........
  // Clean up
  // ........
  delete[] y; delete[] dydp; 
}

// ###################################
// Right hand sides of layer equations
// ###################################
void FourField::CashKarp45Rhs (double x, complex<double>* y, complex<double>* dydx)
{
  if (rhs_chooser == 0)
    {
      // ...............................
      // Define layer equation variables
      // ...............................
      complex<double> g (g_r, g_i);
      
      double p  = x;
      double p2 = p*p;
      double p3 = p*p2;
      double p4 = p2*p2;
      double D2 = D*D;

      // ........................................................................
      // Right-hand sides for backward integration of three-field layer equations
      // ........................................................................
      complex<double> W   = y[0];
      complex<double> V   = y[1];
      complex<double> gEe = g + mc_i * Qe;
      complex<double> gEi = g + mc_i * Qi;
      complex<double> gPD = Pperp + (g + mc_i * Qi) * D*D;
      double          PS  = Pphi + Pperp;
      double          PP  = Pphi * Pperp;
      double          PD  = Pphi * D*D /iotae;
      
      complex<double> A = p2 /(gEe + p2); 
      complex<double> B = g * gEi + gEi * PS * p2 + PP * p4;
      complex<double> C;
      if (!lowD)
	{
	  C = gEe + gPD * p2 + PD * p4;
	}
      else
	{
	  C = gEe + Pperp * p2;
	}
      
      complex<double> AA = (gEe - p2) /(gEe + p2);
      
      dydx[0] = - AA * W /p - W*W /p + B * p3 /A /C;
    }
  else
    {
      // ...............................
      // Define layer equation variables
      // ...............................
      complex<double> g (g_r, g_i);
      
      double p  = x;
      double p2 = p*p;
      double p4 = p2*p2;
      double D2 = D*D;
      double cm = 1./cbeta/cbeta;

      complex<double> gEe = g + mc_i * Qe;
      complex<double> gEi = g + mc_i * Qi;
      complex<double> gPD = Pperp + (g + mc_i * Qi) * D*D;
      complex<double> gP2 = g + Pphi * p2;

      complex<double> E11 = 2. * gEe /(gEe + p2);
      complex<double> E21 = - 2. * mc_i * Qe * (g + 2. * Pphi * p2) /(gEe + p2) /gP2;
      complex<double> E22 = - 2. * Pphi * p2 /gP2;

      complex<double> F11 = p2 * (gEe + p2) * (gEi + Pphi * p2);
      complex<double> F12 = p2 * (gEe + p2) * (gEi + Pphi * p2 /iotae);
      complex<double> F21 = - mc_i * Qe * p2 * (gEi + Pphi * p2)
	+ cm * p2 * gP2 * (mc_i * Qe + D2 * gEi * p2 + D2 * Pphi * p4);
      complex<double> F22 = - mc_i * Qe * p2 * (gEi + Pphi * p2 /iotae)
	+ cm * p2 * gP2 * (gEe + gPD * p2 + D2 * Pphi * p4 /iotae);

      complex<double> W11 = y[0];
      complex<double> W12 = y[1];
      complex<double> W21 = y[2];
      complex<double> W22 = y[3];

      // .......................................................................
      // Right-hand sides for backward integration of four-field layer equations
      // .......................................................................
      complex<double> rhs11 = W11 - W11 * W11 - W12 * W21 - E11 * W11             + F11;
      complex<double> rhs12 = W12 - W11 * W12 - W12 * W22 - E11 * W12             + F12;
      complex<double> rhs21 = W21 - W21 * W11 - W22 * W21 - E21 * W11 - E22 * W21 + F21;
      complex<double> rhs22 = W22 - W21 * W12 - W22 * W22 - E21 * W12 - E22 * W22 + F22;

      dydx[0] = rhs11 /p;
      dydx[1] = rhs12 /p;
      dydx[2] = rhs21 /p;
      dydx[3] = rhs22 /p;
    }
}

// #########################################################
// Function to calculate Jace Waybright's approximate Deltas
// #########################################################
void FourField::GetJace (complex<double>& _DeltaJ1, complex<double>& _DeltaJ2)
{
  complex<double> g (g_r, g_i);
  complex<double> ge = g + complex<double> (0., Qe);
  complex<double> gi = g + complex<double> (0., Qi);

  double G14 = gsl_sf_gamma (0.25);
  double G34 = gsl_sf_gamma (0.75);
  double G13 = gsl_sf_gamma (1./3.);
  double G53 = gsl_sf_gamma (5./3.);
  double G56 = gsl_sf_gamma (5./6.);
  double G76 = gsl_sf_gamma (7./6.);
  
  complex<double> a1 = - cbeta*cbeta /g /ge;
  complex<double> a2 = - (2. * G34 /G14) * complex<double> (- 1./sqrt(2.), 1./sqrt(2.)) * pow (g * gi /ge, 0.25);
  complex<double> a3 = - a2 * (ge /a1 + 6. * a1 * g * gi);
  complex<double> a4 = 1. - 2. * sqrt (a1) * g * gi /a3;
  complex<double> a5 = 1. /a2 + 6. * a1 * g * gi /a3;

  _DeltaJ1 = ge * a4 /a5;

  complex<double> b1 = g * g * gi /cbeta/cbeta;
  complex<double> b2 = pow (b1,  1./24.) * G13 /pow (2., 5./3.) /pow (3., 1./6.) /sqrt (M_PI);
  complex<double> b3 = pow (b1, -1./8.)  * sqrt (3.) * G13 * G76 /4. /M_PI;
  complex<double> b4 = pow (b1, -5./8.)  * 9. * sqrt (3.) * G56 * G53 /2. /M_PI;
  complex<double> b5 = b2 - b4 * g * gi /12. /ge;

  _DeltaJ2 = - M_PI * ge * b5 /b3;
}
