#include "Layer.h"

// ###########
// Constructor
// ###########
Layer::Layer ()
{
  // ---------------------------------------------
  // Ensure that directory ../Outputs/Layer exists
  // ---------------------------------------------
  if (!CreateDirectory ("../Outputs"))
    {
      exit (1);
    }
  if (!CreateDirectory ("../Outputs/Layer"))
    {
      exit (1);
    }
  
  // ............................
  // Set miscellaneous parameters
  // ............................
  Im = complex<double> (0., 1.);

  // .........................................
  // Read control parameters from TJ JSON file
  // .........................................
  string JSONFilename = "../Inputs/TJ.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  RMP   = JSONData["RMP"]  .get<int> ();

  // ......................................
  // Read control parameters from JSON file
  // ......................................
  string JSONFilename1 = "../Inputs/Layer.json";
  json   JSONData1     = ReadJSONFile (JSONFilename1);

  MARG    = JSONData1["MARG"]  .get<int>    ();

  pstart  = JSONData1["pstart"].get<double> ();
  pend    = JSONData1["pend"]  .get<double> ();
  P3max   = JSONData1["P3max"] .get<double> ();
  Nscan   = JSONData1["Nscan"] .get<int>    ();

  acc  = JSONData1["acc"] .get<double> ();
  h0   = JSONData1["h0"]  .get<double> ();
  hmin = JSONData1["hmin"].get<double> ();
  hmax = JSONData1["hmax"].get<double> ();

  dS      = JSONData1["dS"]     .get<double> ();
  Smax    = JSONData1["Smax"]   .get<double> ();
  Smin    = JSONData1["Smin"]   .get<double> ();
  Eps     = JSONData1["Eps"]    .get<double> ();
  MaxIter = JSONData1["MaxIter"].get<int>    ();

  // ------------
  // Sanity check
  // ------------
  if (pend < 0.)
    {
      printf ("Layer:: Error - pstart cannot be negative\n");
      exit (1);
    }
  if (pstart < pend)
    {
      printf ("Layer:: Error - pstart cannot be less than pend\n");
      exit (1);
    }
  if (P3max < 0.)
    {
      printf ("Layer:: Error - P3max cannot be negative\n");
      exit (1);
    }
  if (Nscan < 0)
    {
      printf ("Layer:: Error - Nscan cannot be negative\n");
      exit (1);
    }
  if (acc < 0.)
    {
      printf ("Layer:: Error - acc cannot be negative\n");
      exit (1);
    }
  if (h0 < 0.)
    {
      printf ("Layer:: Error - h0 cannot be negative\n");
      exit (1);
    }
   if (hmin < 0.)
    {
      printf ("Layer:: Error - hmin cannot be negative\n");
      exit (1);
    }
   if (hmax < hmin)
    {
      printf ("Layer:: Error - hmax cannot be less than hmin\n");
      exit (1);
    }
   if (dS < 0.)
    {
      printf ("Layer:: Error - dS cannot be negative\n");
      exit (1);
    }
  if (Smin < 0.)
    {
      printf ("Layer:: Error - Smin cannot be negative\n");
      exit (1);
    }
  if (Smax < Smin)
    {
      printf ("Layer:: Error - Smax cannot be less than Smin\n");
      exit (1);
    }
  if (Eps < 0.)
    {
      printf ("Layer:: Error - Eps cannot be less than zero\n");
      exit (1);
    }
   if (Maxiter < 0)
    {
      printf ("Layer:: Error - Maxiter cannot be less than zero\n");
      exit (1);
    }
    
  printf ("\n");
  printf ("Class LAYER::\n");
  printf ("Git Hash     = "); printf (GIT_HASH);     printf ("\n");
  printf ("Compile time = "); printf (COMPILE_TIME); printf ("\n");
  printf ("Git Branch   = "); printf (GIT_BRANCH);   printf ("\n\n");
  printf ("pstart = %10.3e pend = %10.3e P3max = %10.3e Nscan =  %-4d      MARG    = %-1d\n",
	  pstart, pend, P3max, Nscan, MARG);
  printf ("acc    = %10.3e h0   = %10.3e hmin  = %10.3e hmax  = %10.3e\n",
	  acc, h0, hmin, hmax);
  printf ("dS     = %10.3e Smax = %10.3e Smin  = %10.3e Eps   = %10.3e MaxIter = %-4d\n",
	  dS, Smax, Smin, Eps, MaxIter);
}

// #########################
// Function to solve problem
// #########################
void Layer::Solve (int verbose)
{
  // .....................................
  // Read in rational surface data from TJ
  // .....................................
  ReadNetcdf ();

  // ...............
  // Allocate memory
  // ...............
  np_marg.resize (nres);
  gr_marg.resize (nres, 10);
  gi_marg.resize (nres, 10);
  Dr_marg.resize (nres, 10);
  Di_marg.resize (nres, 10);

  for (int i = 0; i < nres; i++)
    for (int j = 0; j < 10; j++)
      {
	gr_marg(i, j) = - 1.e15;
	gi_marg(i, j) = - 1.e15;
	Dr_marg(i, j) = - 1.e15;
	Di_marg(i, j) = - 1.e15;
      }

  gamma_e = new double[nres];
  omega_e = new double[nres];
  f_e     = new double[nres];
  res_e   = new double[nres];
  lowD_e  = new int   [nres];

  omega_r.resize (nres, Nscan + 1);
  Deltar.resize  (nres, Nscan + 1);
  Deltai.resize  (nres, Nscan + 1);
  Xi_res.resize  (nres, Nscan + 1);
  T_res.resize   (nres, Nscan + 1);

  // ..............................
  // Find marginal stability points
  // ..............................
  printf ("Marginal stability points:\n");
  for (int i = 0; i < nres; i++)
    FindMarginal (i, verbose);
  
  // ...........................................................
  // Calculate electron-branch growth-rates and real frequencies
  // ...........................................................
  printf ("Electron-branch growth-rates and real frequencies:\n");
  for (int i = 0; i < nres; i++)
    GetElectronBranchGrowth (i, verbose);
  
  // ............................................
  // Calculate shielding factor and torque curves
  // ............................................
  if (RMP)
    {
      printf ("Calculating shielding factor and torque curves:\n");
      for (int i = 0; i < nres; i++)
	GetTorque (i);
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

// #########################################
// Function to read TJ data from netcdf file
// #########################################
void Layer::ReadNetcdf ()
{
  printf ("Reading data from netcdf file Outputs/TJ/TJ.nc:\n");

  try
    {
      NcFile dataFile ("../Outputs/TJ/TJ.nc", NcFile::read);

      NcVar input_x  = dataFile.getVar ("InputParameters");
      NcVar rres_x   = dataFile.getVar ("r_res");
      NcVar mpol_x   = dataFile.getVar ("m_res");
      NcVar S13_x    = dataFile.getVar ("S13");
      NcVar tau_x    = dataFile.getVar ("tau");
      NcVar QE_x     = dataFile.getVar ("QE");
      NcVar Qe_x     = dataFile.getVar ("Qe");
      NcVar Qi_x     = dataFile.getVar ("Qi");
      NcVar iotae_x  = dataFile.getVar ("iota_e");
      NcVar D_x      = dataFile.getVar ("D");
      NcVar Pphi_x   = dataFile.getVar ("Pphi");
      NcVar Pperp_x  = dataFile.getVar ("Pperp");
      NcVar Delta_x  = dataFile.getVar ("Delta");
      NcVar Deltac_x = dataFile.getVar ("Delta_crit");
      NcDim n_x      = rres_x.getDim (0);
      NcDim p_x      = input_x.getDim (0);

      NcVar Chia_x;
      if (RMP)
	  Chia_x = dataFile.getVar ("Chi_a");
 
      nres       = n_x.getSize ();
      int npara  = p_x.getSize ();
      input      = new double[npara];
      r_res      = new double[nres];
      m_res      = new int[nres];
      S13_res    = new double[nres];
      tau_res    = new double[nres];
      QE_res     = new double[nres];
      Qe_res     = new double[nres];
      Qi_res     = new double[nres];
      iotae_res  = new double[nres];
      D_res      = new double[nres];
      Pphi_res   = new double[nres];
      Pperp_res  = new double[nres];
      Delta_res  = new double[nres];
      Deltac_res = new double[nres];
      Chi_res    = new double[nres];

      input_x.getVar  (input);
      rres_x.getVar   (r_res);
      mpol_x.getVar   (m_res);
      S13_x.getVar    (S13_res);
      tau_x.getVar    (tau_res);
      QE_x.getVar     (QE_res);
      Qe_x.getVar     (Qe_res);
      Qi_x.getVar     (Qi_res);
      iotae_x.getVar  (iotae_res);
      D_x.getVar      (D_res);
      Pphi_x.getVar   (Pphi_res);
      Pperp_x.getVar  (Pperp_res);
      Delta_x.getVar  (Delta_res);
      Deltac_x.getVar (Deltac_res);

      if (RMP)
	Chia_x.getVar (Chi_res);
    }
  catch (NcException& e)
     {
       printf ("Error reading data from netcdf file Outputs/TJ/TJ.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }

  printf ("Rational surfaces:\n");
  for (int i = 0; i < nres; i++)
    {
      printf ("m = %3d r = %9.2e De = %9.2e Dc = %9.2e S13 = %9.2e tau = %9.2e Qe = %9.2e Qi = %9.2e D = %9.2e P = %9.2e\n",
	      m_res[i], r_res[i], Delta_res[i], Deltac_res[i], S13_res[i], tau_res[i], Qe_res[i], Qi_res[i],
	      D_res[i], Pphi_res[i]);
    }
}

// ###############################################################################
// Function to find marginal stability points associated with ith rational surface
// ###############################################################################
void Layer::FindMarginal (int i, int verbose)
{
  // ....................
  // Set layer parameters
  // ....................
  Qe    = Qe_res[i];
  Qi    = Qi_res[i];
  D     = D_res[i];
  Pphi  = Pphi_res[i];
  Pperp = Pperp_res[i];
  iotae = iotae_res[i];

  // ---------------------------------------------
  // Set electron branch marginal stability points
  // ---------------------------------------------
  g_r           = 0.;
  np_marg(i)    = 1;
  gr_marg(i, 0) = g_r;
  gi_marg(i, 0) = - Qe;
  Dr_marg(i, 0) = g_r;
  Di_marg(i, 0) = g_r;

  // .........................................
  // Find ion branch marginal stability points
  // .........................................
  if (MARG)
    {
      // ..........................................
      // Perform frequency scan at zero growth-rate
      // ..........................................
      double* gg_i = new double[Nscan + 1];
      double* DD_i = new double[Nscan + 1];
      
      for (int j = 1; j <= Nscan; j++)
	{
	  g_i = - QE - Qe - (Qi - Qe) * double (j) /double (Nscan);
	  
	  SolveLayerEquations ();
	  gg_i[j] = g_i;
	  DD_i[j] = imag (Deltas);
	  
	  if (verbose & (j-1)%100 == 0)
	    printf ("Surface %1d: g_r = %11.4e g_i = %11.4e Deltas_i = %11.4e\n",
		    i+1, g_r, g_i, DD_i[j]);
	}
       
      int cnt = 0;
      for (int j = 1; j < Nscan; j++)
	{
	  double x;
	  if (DD_i[j] * DD_i[j+1] < 0.)
	    {
	      Ridder (gg_i[j], gg_i[j+1], DD_i[j], DD_i[j+1], x);
	      
	      g_i = x;
	      SolveLayerEquations ();
	      
	      cnt++;
	      np_marg(i)     += 1;
	      gr_marg(i, cnt) = g_r;
	      gi_marg(i, cnt) = g_i;
	      Dr_marg(i, cnt) = real (Deltas);
	      Di_marg(i, cnt) = imag (Deltas);
	    }
	}
      delete[] gg_i; delete[] DD_i;
    }

  for (int j = 0; j < np_marg(i); j++)
    printf ("Rational surface %3d: g = (%10.3e, %10.3e) Delta_s = (%10.3e, %10.3e)\n",
	    i+1, gr_marg(i, j), gi_marg(i, j), Dr_marg(i, j), Di_marg(i, j));
}

// #################################################################################################
// Function to find complex growth-rate of electron branch-mode associated with ith rational surface
// #################################################################################################
void Layer::GetElectronBranchGrowth (int i, int verbose)
{
  // ....................
  // Set layer parameters
  // ....................
  Qe    = Qe_res[i];
  Qi    = Qi_res[i];
  QE    = QE_res[i];
  D     = D_res[i];
  Pphi  = Pphi_res[i];
  Pperp = Pperp_res[i];
  iotae = iotae_res[i];

  // ..........................................
  // Search for solution of dispersion relation
  // ..........................................
  Delta = (Delta_res[i] - Deltac_res[i]) /S13_res[i];

  double gr, gi, Residual;
  gr = 0.;
  gi = gi_marg (i, 0);
  NewtonRoot (gr, gi, Residual, verbose);

  double f     =    (gi) * (gi + Qi) /(- Qe) /(- Qe + Qi);
  f            += - (gi) * (gi + Qe) /(- Qi) /(- Qi + Qe);
  double gamma = gr /tau_res[i];
  double omega = (QE - g_i) /tau_res[i];

  gamma_e[i] = gamma/1.e3;
  omega_e[i] = omega/1.e3;
  f_e    [i] = f;
  res_e  [i] = Residual;
  lowD_e [i] = lowD;
  
  printf ("Rational surface %3d: gamma = %10.3e (kHz) omega = %10.3e (kHz) f = %10.3e res = %10.3e lowD = %1d\n",
	  i+1, gamma/1.e3, omega/1.e3, f, Residual, lowD);
}

// #############################################################################################
// Function to calculate shielding factor and torque curves associated with ith rational surface
// #############################################################################################
void Layer::GetTorque (int i)
{
  // ....................
  // Set layer parameters
  // ....................
  Qe    = Qe_res[i];
  Qi    = Qi_res[i];
  D     = D_res[i];
  Pphi  = Pphi_res[i];
  Pperp = Pperp_res[i];
  iotae = iotae_res[i];

  // .......................................
  // Calculate torques and shielding factors
  // .......................................
  double Dr, Di, A, B;
  
  for (int j = 0; j <= Nscan; j++)
    omega_r(i ,j) = (1.1*Qi_res[i] + QE_res[i] + (2.*Qe_res[i] - 1.1*Qi_res[i] + QE_res[i]) * double (j) /double (Nscan)) /tau_res[i];
  
  for (int j = 0; j <= Nscan; j++)
    {
      g_r = 0.;
      g_i = QE_res[i] - omega_r(i, j) * tau_res[i];

      SolveLayerEquations ();

      B  = Deltac_res[i] - Delta_res[i] ;
      Dr = S13_res[i] * real (Deltas) + B;
      Di = S13_res[i] * imag (Deltas);
      A  = sqrt (Dr*Dr + Di*Di);
      
      Deltar(i, j) = Dr;
      Deltai(i, j) = Di;

      Xi_res(i, j) = fabs (B) /A;

      T_res(i, j)  = 2.*M_PI*M_PI * input[0] * Chi_res[i]*Chi_res[i] * Di /A/A;
    }
}

// ###########################################
// Function to write Layer data to netcdf file
// ###########################################
void Layer::WriteNetcdf ()
{
   printf ("Writing data to netcdf file Outputs/Layer/Layer.nc:\n");

   int*    np_y = new int[nres];
   double* gr_y = new double[nres*10];
   double* gi_y = new double[nres*10];
   double* dr_y = new double[nres*10];
   double* di_y = new double[nres*10];

   double* om_y = new double[nres*(Nscan+1)];
   double* Dr_y = new double[nres*(Nscan+1)];
   double* Di_y = new double[nres*(Nscan+1)];
   double* xi_y = new double[nres*(Nscan+1)];
   double* t_y  = new double[nres*(Nscan+1)];

   for (int i = 0; i < nres; i++)
     np_y[i] = np_marg(i);

   int cnt = 0;
   for (int i = 0; i < nres; i++)
     for (int j = 0; j < 10; j++)
       {
	 gr_y[cnt] = gr_marg(i, j);
	 gi_y[cnt] = gi_marg(i, j);
	 dr_y[cnt] = Dr_marg(i, j);
	 di_y[cnt] = Di_marg(i, j);
	 cnt++;
       }

   cnt = 0;
   for (int i = 0; i < nres; i++)
     for (int j = 0; j <= Nscan; j++)
       {
	 om_y[cnt] = omega_r(i, j);
	 Dr_y[cnt] = Deltar(i, j);
	 Di_y[cnt] = Deltai(i, j);
	 xi_y[cnt] = Xi_res(i, j);
	 t_y[cnt]  = T_res(i, j);
	 cnt++;
       }
      
   try
     {
       NcFile dataFile ("../Outputs/Layer/Layer.nc", NcFile::replace);

       dataFile.putAtt ("Git_Hash",     GIT_HASH);
       dataFile.putAtt ("Compile_Time", COMPILE_TIME);
       dataFile.putAtt ("Git_Branch",   GIT_BRANCH);

       NcDim x_d = dataFile.addDim ("nres",  nres);
       NcDim y_d = dataFile.addDim ("nmarg", 10);
       NcDim z_d = dataFile.addDim ("Nscan", Nscan+1);

       vector<NcDim> marg_d;
       marg_d.push_back (x_d);
       marg_d.push_back (y_d);

       vector<NcDim> torq_d;
       torq_d.push_back (x_d);
       torq_d.push_back (z_d);

       NcVar rres_x      = dataFile.addVar ("r_res",      ncDouble, x_d);
       rres_x.putVar (r_res);
       NcVar mres_x      = dataFile.addVar ("m_res",      ncInt,    x_d);
       mres_x.putVar (m_res);
       NcVar Deltares_x  = dataFile.addVar ("Delta_res",  ncDouble, x_d);
       Deltares_x.putVar (Delta_res);
       NcVar Deltacres_x = dataFile.addVar ("Deltac_res", ncDouble, x_d);
       Deltacres_x.putVar (Deltac_res);
       NcVar S13res_x    = dataFile.addVar ("S13_res",    ncDouble, x_d);
       S13res_x.putVar (S13_res);
       NcVar taures_x    = dataFile.addVar ("tau_res",    ncDouble, x_d);
       taures_x.putVar (tau_res);
       NcVar QEres_x     = dataFile.addVar ("QE_res",     ncDouble, x_d);
       QEres_x.putVar (QE_res);
       NcVar Qeres_x     = dataFile.addVar ("Qe_res",     ncDouble, x_d);
       Qeres_x.putVar (Qe_res);
       NcVar Qires_x     = dataFile.addVar ("Qi_res",     ncDouble, x_d);
       Qires_x.putVar (Qi_res);
       NcVar iotaeres_x  = dataFile.addVar ("iotae_res",  ncDouble, x_d);
       iotaeres_x.putVar (iotae_res);
       NcVar D_x         = dataFile.addVar ("D_res",      ncDouble, x_d);
       D_x.putVar (D_res);
       NcVar Pphires_x   = dataFile.addVar ("Pphi_res",   ncDouble, x_d);
       Pphires_x.putVar (Pphi_res);
       NcVar Pperpres_x  = dataFile.addVar ("Pperp_res",  ncDouble, x_d);
       Pperpres_x.putVar (Pperp_res);
       
       NcVar gammae_x = dataFile.addVar ("gamma_e", ncDouble, x_d);
       gammae_x.putVar (gamma_e);
       NcVar omegae_x = dataFile.addVar ("omega_e", ncDouble, x_d);
       omegae_x.putVar (omega_e);
       NcVar fe_x     = dataFile.addVar ("f_e",     ncDouble, x_d);
       fe_x.putVar (f_e);
       NcVar rese_x   = dataFile.addVar ("res_e",   ncDouble, x_d);
       rese_x.putVar (res_e);
       NcVar lowDe_x  = dataFile.addVar ("lowD_e",  ncInt,    x_d);
       lowDe_x.putVar (lowD_e);

       if (MARG)
	 {
	   NcVar np_x = dataFile.addVar ("n_marg",  ncInt,    x_d);
	   np_x.putVar (np_y);
	   NcVar gr_x = dataFile.addVar ("gr_marg", ncDouble, marg_d);
	   gr_x.putVar (gr_y);
	   NcVar gi_x = dataFile.addVar ("gi_marg", ncDouble, marg_d);
	   gi_x.putVar (gi_y);
	   NcVar Dr_x = dataFile.addVar ("Dr_marg", ncDouble, marg_d);
	   Dr_x.putVar (dr_y);
	   NcVar Di_x = dataFile.addVar ("Di_marg", ncDouble, marg_d);
	   Di_x.putVar (di_y);
	 }

       if (RMP)
	 {
	   NcVar om_x = dataFile.addVar  ("omega_r", ncDouble, torq_d);
	   om_x.putVar (om_y);
	   NcVar DDr_x = dataFile.addVar ("Deltar",  ncDouble, torq_d);
	   DDr_x.putVar (Dr_y);
	   NcVar DDi_x = dataFile.addVar ("Deltai",  ncDouble, torq_d);
	   DDi_x.putVar (Di_y);
	   NcVar xi_x = dataFile.addVar  ("Xi_res",  ncDouble, torq_d);
	   xi_x.putVar (xi_y);
	   NcVar t_x  = dataFile.addVar  ("T_res",   ncDouble, torq_d);
	   t_x.putVar (t_y);
	 }
     }
   catch (NcException& e)
     {
       printf ("Error writing data to netcdf file Outputs/Layer/Layer.nc\n");
       printf ("%s\n", e.what ());
       exit (1);
     }

   delete[] gr_y; delete[] gi_y; delete[] dr_y; delete[] di_y;
   delete[] om_y; delete[] xi_y; delete[] t_y;
}

// #############################
// Function to deallocate memory
// #############################
void Layer::CleanUp ()
{
  delete[] r_res,    delete[] m_res;     delete[] Delta_res; delete[] Deltac_res; delete[] tau_res;
  delete[] QE_res;   delete[] Qe_res;    delete[] Qi_res;    delete[] iotae_res;  delete[] D_res;
  delete[] Pphi_res; delete[] Pperp_res; delete[] S13_res;   delete[] Chi_res;    delete[] input;
  delete[] gamma_e;  delete[] omega_e;   delete[] res_e;     delete[] lowD_e;     delete[] f_e;
}

// #################################
// Function to solve layer equations
// #################################
void Layer::SolveLayerEquations ()
{
  // ........
  // Define g
  // ........
  complex<double> g = complex<double> (g_r, g_i);

  // ..............
  // Determine Pmax
  // ..............
  complex<double> gEe = g + Im * Qe;
  complex<double> gEi = g + Im * Qi;
  complex<double> gPD = Pperp + (g + Im * Qi) * D*D;
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
      lowD = 1;
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
  int               rept; count = 0;
  complex<double>*  y    = new complex<double>[1];
  complex<double>*  dydp = new complex<double>[1];

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
      gamma = - Im * (Qe - Qi) * Pphi /Pperp + gEi;
      X     = (alpha * beta - gamma) /2./ sqrt(beta);

      y[0] = - 1. + X * p - sqrt(beta) * p*p*p;
    }
  
  do
    {
      CashKarp45Adaptive (1, p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
    }
  while (p > pend);
  
  CashKarp45Rhs (p, y, dydp);
  Deltas = M_PI /dydp[0];

  // ........
  // Clean up
  // ........
  delete[] y; delete[] dydp; 
}

// ###################################
// Right hand sides of layer equations
// ###################################
void Layer::CashKarp45Rhs (double x, complex<double>* y, complex<double>* dydx)
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
 
  // ............................................................
  // Right-hand sides for backward integration of layer equations
  // ............................................................
  complex<double> W   = y[0];
  complex<double> V   = y[1];
  complex<double> gEe = g + Im * Qe;
  complex<double> gEi = g + Im * Qi;
  complex<double> gPD = Pperp + (g + Im * Qi) * D*D;
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

// ######################################################################
// Function to calculate target functions for Newton-Raphson root finding
// ######################################################################
void Layer::NewtonFunction (double x1, double x2, double& F1, double& F2)
{
  g_r = x1;
  g_i = x2;

  SolveLayerEquations ();

  F1 = real (Deltas) - Delta;
  F2 = imag (Deltas);
}

// ################################
// Target function for root finding
// ################################
double Layer::RootFindF (double x)
{
  g_i = x;

  SolveLayerEquations ();

  return imag (Deltas);
}

