#include "Layer.h"

// ###########
// Constructor
// ###########
Layer::Layer ()
{
  // ...................................
  // Set adaptive integration parameters
  // ...................................
  maxrept = 50;
  flag    = 2;
  
  // ........................
  // Set Cash-Karp parameters
  // ........................
  aa1 = 0.;
  aa2 = 1. /5.;
  aa3 = 3. /10.;
  aa4 = 3. /5.;
  aa5 = 1.;
  aa6 = 7. /8.;
  
  cc1 = 37.  /378.;
  cc3 = 250. /621.;
  cc4 = 125. /594.;
  cc6 = 512. /1771.;
  
  ca1 = cc1 - 2825.  /27648.;
  ca3 = cc3 - 18575. /48384.;
  ca4 = cc4 - 13525. /55296.;
  ca5 =     - 277.   /14336.;
  ca6 = cc6 - 1.     /4.;
  
  bb21 = 1. /5.;
  
  bb31 = 3. /40.;
  bb32 = 9. /40.;
  
  bb41 =   3. /10.;
  bb42 = - 9. /10.;
  bb43 =   6. /5.;
  
  bb51 = - 11. /54.;
  bb52 =    5. /2.;
  bb53 = - 70. /27.;
  bb54 =   35. /27.;
  
  bb61 = 1631.  /55296.;
  bb62 = 175.   /512.;
  bb63 = 575.   /13824.;
  bb64 = 44275. /110592.;
  bb65 = 253.   /4096.;

  // ............................
  // Set miscellaneous parameters
  // ............................
  Im = complex<double> (0., 1.);

  // ......................................
  // Read control parameters from JSON file
  // ......................................
  string JSONFilename = "Inputs/Layer.json";
  json   JSONData     = ReadJSONFile (JSONFilename);

  LAYER   = JSONData["LAYER"] .get<int>    ();
  pstart  = JSONData["pstart"].get<double> ();
  pend    = JSONData["pend"]  .get<double> ();
  P3max   = JSONData["P3max"] .get<double> ();
  Nscan   = JSONData["Nscan"] .get<int>    ();

  acc  = JSONData["acc"] .get<double> ();
  h0   = JSONData["h0"]  .get<double> ();
  hmin = JSONData["hmin"].get<double> ();
  hmax = JSONData["hmax"].get<double> ();

  eps     = JSONData["eps"]    .get<double> ();
  smax    = JSONData["smax"]   .get<double> ();
  smin    = JSONData["smin"]   .get<double> ();
  Eta     = JSONData["Eta"]    .get<double> ();
  Maxiter = JSONData["Maxiter"].get<int>    ();

  if (LAYER)
    {
      printf ("\n");
      printf ("Class LAYER::\n");
      printf ("pstart = %10.3e pend = %10.3e P3max = %10.3e Nscan = %4d\n",
	      pstart, pend, P3max, Nscan);
      printf ("acc    = %10.3e h0   = %10.3e hmin  = %10.3e hmax  = %10.3e\n",
	      acc, h0, hmin, hmax);
      printf ("eps    = %10.3e smax = %10.3e smin  = %10.3e Eta   = %10.3e Maxiter = %4d\n",
	      eps, smax, smin, Eta, Maxiter);
    }
}

// #########################
// Function to solve problem
// #########################
void Layer::Solve ()
{
  // ....................................
  // Skip calculation unless LAYER is set
  // ....................................
  if (!LAYER)
    return;
  
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
    FindMarginal (i);
  
  // ...........................................................
  // Calculate electron-branch growth-rates and real frequencies
  // ...........................................................
  printf ("Electron-branch growth-rates and real frequencies:\n");
  for (int i = 0; i < nres; i++)
    GetElectronBranchGrowth (i);

  // ............................................
  // Calculate shielding factor and torque curves
  // ............................................
  printf ("Calculating shielding factor and torque curves:\n");
  for (int i = 0; i < nres; i++)
    GetTorque (i);
     
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
      NcFile dataFile ("Outputs/TJ/TJ.nc", NcFile::read);

      NcVar input_x  = dataFile.getVar ("InputParameters");
      NcVar rres_x   = dataFile.getVar ("r_res");
      NcVar mpol_x   = dataFile.getVar ("m_res");
      NcVar Delta_x  = dataFile.getVar ("Delta");
      NcVar Deltac_x = dataFile.getVar ("Delta_crit");
      NcVar Chia_x   = dataFile.getVar ("Chi_a");
      NcVar S13_x    = dataFile.getVar ("S13");
      NcVar tau_x    = dataFile.getVar ("tau");
      NcVar QE_x     = dataFile.getVar ("QE");
      NcVar Qe_x     = dataFile.getVar ("Qe");
      NcVar Qi_x     = dataFile.getVar ("Qi");
      NcVar iotae_x  = dataFile.getVar ("iota_e");
      NcVar D_x      = dataFile.getVar ("D");
      NcVar Pphi_x   = dataFile.getVar ("Pphi");
      NcVar Pperp_x  = dataFile.getVar ("Pperp");
      NcDim n_x      = rres_x.getDim (0);
      NcDim p_x      = input_x.getDim (0);

      nres       = n_x.getSize ();
      int npara  = p_x.getSize ();
      input      = new double[npara];
      r_res      = new double[nres];
      m_res      = new int[nres];
      Delta_res  = new double[nres];
      Deltac_res = new double[nres];
      Chi_res    = new double[nres];
      S13_res    = new double[nres];
      tau_res    = new double[nres];
      QE_res     = new double[nres];
      Qe_res     = new double[nres];
      Qi_res     = new double[nres];
      iotae_res  = new double[nres];
      D_res      = new double[nres];
      Pphi_res   = new double[nres];
      Pperp_res  = new double[nres];

      input_x.getVar  (input);
      rres_x.getVar   (r_res);
      mpol_x.getVar   (m_res);
      Delta_x.getVar  (Delta_res);
      Deltac_x.getVar (Deltac_res);
      Chia_x.getVar   (Chi_res);
      S13_x.getVar    (S13_res);
      tau_x.getVar    (tau_res);
      QE_x.getVar     (QE_res);
      Qe_x.getVar     (Qe_res);
      Qi_x.getVar     (Qi_res);
      iotae_x.getVar  (iotae_res);
      D_x.getVar      (D_res);
      Pphi_x.getVar   (Pphi_res);
      Pperp_x.getVar  (Pperp_res);
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
void Layer::FindMarginal (int i)
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

  // ..........................................
  // Perform frequency scan at zero growth-rate
  // ..........................................
  double* gg_i = new double[Nscan + 1];
  double* DD_i = new double[Nscan + 1];
  
  g_r = 0.;
  for (int j = 1; j <= Nscan; j++)
    {
      g_i = - Qe - (Qi - Qe) * double (j) /double (Nscan);

      SolveLayerEquations ();
      gg_i[j] = g_i;
      DD_i[j] = imag (Deltas);
    }

  // ..............................
  // Find marginal stability points
  // ..............................
  int cnt = 0;
  np_marg(i)    = 1;
  gr_marg(i, 0) = g_r;
  gi_marg(i, 0) = - Qe;
  Dr_marg(i, 0) = g_r;
  Di_marg(i, 0) = g_r;

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

  for (int j = 0; j < np_marg(i); j++)
    printf ("Rational surface %3d: g = (%10.3e, %10.3e) Delta_s = (%10.3e, %10.3e)\n",
	    i+1, gr_marg(i, j), gi_marg(i, j), Dr_marg(i, j), Di_marg(i, j));
}

// #################################################################################################
// Function to find complex growth-rate of electron branch-mode associated with ith rational surface
// #################################################################################################
void Layer::GetElectronBranchGrowth (int i)
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

  // ..........................................
  // Search for solution of dispersion relation
  // ..........................................
  Delta = (Delta_res[i] - Deltac_res[i]) /S13_res[i];

  double gr, gi, F;
  gr = 0.;
  gi = gi_marg (i, 0);
  GetRoot (gr, gi, F);
  
  double gamma = gr /tau_res[i];
  double omega = (QE_res[i] - g_i) /tau_res[i];

  gamma_e[i] = gamma/1.e3;
  omega_e[i] = omega/1.e3;
  res_e  [i] = F;
  lowD_e [i] = lowD;
  
  printf ("Rational surface %3d: gamma = %10.3e (kHz) omega = %10.3e (kHz) res = %10.3e lowD = %1d\n",
	  i+1, gamma/1.e3, omega/1.e3, F, lowD);
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
       NcFile dataFile ("Outputs/Layer/Layer.nc", NcFile::replace);

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
       NcVar rese_x   = dataFile.addVar ("res_e",   ncDouble, x_d);
       rese_x.putVar (res_e);
       NcVar lowDe_x  = dataFile.addVar ("lowD_e",  ncInt,    x_d);
       lowDe_x.putVar (lowD_e);

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
  delete[] gamma_e;  delete[] omega_e;   delete[] res_e;     delete[] lowD_e;
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
      CashKarp45Adaptive (1, p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag);
    }
  while (p > pend);
  
  Rhs (p, y, dydp);
  Deltas = M_PI /dydp[0];

  // ........
  // Clean up
  // ........
  delete[] y; delete[] dydp; 
}

// ###################################
// Right hand sides of layer equations
// ###################################
void Layer::Rhs (double x, complex<double>* y, complex<double>* dydx)
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

// #######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step-length Cash-Karp fourth-order/fifth-order
//  Runge-Kutta scheme
//
//     neqns   ... number of coupled equations
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step-length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... safety factor
//     T       ... step-length cannot change by more than this factor from step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step-length
//     h_max   ... maximum allowable step-length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step while attempting to maintain 
//  constant truncation error per step of acc:
//
//     flag = 0 ... error is absolute
//     flag = 1 ... error is relative
//     flag = 2 ... error is mixed
//
// #######################################################################
void Layer::CashKarp45Adaptive (int neqns, double& x, complex<double>* y, double& h, 
				double& t_err, double acc, double S, double T, int& rept,
				int maxrept, double h_min, double h_max, int flag)
{
  complex<double>* y0  = new complex<double>[neqns];
  complex<double>* Err = new complex<double>[neqns];
  double           hin = h;

  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0[i] = y[i];

  // Take Cash-Karp RK4/RK5 step 
  CashKarp45Fixed (neqns, x, y, Err, h);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = abs (Err[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = abs (Err[i] /y[i]);
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = abs (Err[i] /y[i]);
	  err2  = abs (Err[i]);
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step-length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc /t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc /t_err), 0.25);

  // Prevent step-length from changing by more than factor T
  if (h_est /h > T)
    h *= T;
  else if (h_est /h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h /fabs(h) : h;

  // Prevent step-length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = + h_min;
      else
	h = - h_min;
    }

  // Check if truncation error acceptable
  if ((t_err <= acc) || (count >= maxrept))
    {
      // If truncation error acceptable take step 
      rept  = count;
      count = 0;
    }
  else 
    {
      // If truncation error unacceptable repeat step 
      count++;
      x = x0;
      for (int i = 0; i < neqns; i++)
	y[i] = y0[i];
      CashKarp45Adaptive (neqns, x, y, h, t_err, acc, S, T, rept, 
			  maxrept, h_min, h_max, flag);
    }

  delete[] y0; delete[] Err;
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step-length Cash-Karp fourth-order/fifth-order
// Runge-Kutta scheme
//
//     neqns ... number of coupled equations
//     x     ... independent variable
//     y     ... array of dependent variables 
//     err   ... array of errors
//     h     ... step-length
//     
// #####################################################################
void Layer::CashKarp45Fixed (int neqns, double& x, complex<double>* y, complex<double>* err, double h)
{
  complex<double>* dydx = new complex<double>[neqns];
  complex<double>* k1   = new complex<double>[neqns];
  complex<double>* k2   = new complex<double>[neqns];
  complex<double>* k3   = new complex<double>[neqns];
  complex<double>* k4   = new complex<double>[neqns];
  complex<double>* k5   = new complex<double>[neqns];
  complex<double>* k6   = new complex<double>[neqns];
  complex<double>* f    = new complex<double>[neqns];

  // First stage
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f [i] = y[i] + bb21 * k1[i];
    }

  // Second stage
  Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f [i] = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }

  // Third stage
  Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f [i] = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }

  // Fourth stage
  Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f [i] = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }

  // Fifth stage
  Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f [i] = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }

  // Sixth stage
  Rhs (x + aa6 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k6[i] = h * dydx[i];
    }

  // Actual step 
  for (int i = 0; i < neqns; i++)
    {
      y  [i] = y[i] + cc1 * k1[i] + cc3 * k3[i] + cc4 * k4[i]               + cc6 * k6[i];
      err[i] =        ca1 * k1[i] + ca3 * k3[i] + ca4 * k4[i] + ca5 * k5[i] + ca6 * k6[i];
    }
  x += h;

  delete[] dydx; delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6;
  delete[] f;
}

// #####################################
// Function to calculate Jacobian matrix
// #####################################
void Layer::GetJacobian (double& J11, double& J12, double& J21, double& J22)
{
  double g_r_save = g_r;
  double g_i_save = g_i;

  g_r = g_r_save - eps;
  g_i = g_i_save;

  SolveLayerEquations ();

  double Deltar1 = real (Deltas);
  double Deltai1 = imag (Deltas);

  g_r = g_r_save + eps;
  g_i = g_i_save;

  SolveLayerEquations ();

  double Deltar2 = real (Deltas);
  double Deltai2 = imag (Deltas);

  J11 = (Deltar2 - Deltar1) /2./eps;
  J21 = (Deltai2 - Deltai1) /2./eps;

  g_r = g_r_save;
  g_i = g_i_save - eps;

  SolveLayerEquations ();

  Deltar1 = real (Deltas);
  Deltai1 = imag (Deltas);

  g_r = g_r_save;
  g_i = g_i_save + eps;

  SolveLayerEquations ();

  Deltar2 = real (Deltas);
  Deltai2 = imag (Deltas);

  J12 = (Deltar2 - Deltar1) /2./eps;
  J22 = (Deltai2 - Deltai1) /2./eps;

  g_r = g_r_save;
  g_i = g_i_save;
}

// ############################################################
// Function to find root of Deltas = Delta via Newton iteration
// ############################################################
 void Layer::GetRoot (double& gr, double& gi, double& F)
{
  double F1, F2, J11, J12, J21, J22, det, iJ11, iJ12, iJ21, iJ22, dgr, dgi, dg, Fold, lambda;
  int    iter, iter1;
  
  g_r = gr;
  g_i = gi;

  SolveLayerEquations ();
  
  F1 = real (Deltas) - Delta;
  F2 = imag (Deltas);

  Fold   = sqrt (F1*F1 + F2*F2);
  lambda = 2.;

  iter = 0;
  do
    {
      iter1 = 0;
      do
	{
	  lambda /= 2.;
	  
	  GetJacobian (J11, J12, J21, J22);
      
	  det = J11 * J22 - J12 * J21;
      
	  iJ11 =   J22 /det;
	  iJ12 = - J12 /det;
	  iJ21 = - J21 /det;
	  iJ22 =   J11 /det;
	  
	  dgr = - (iJ11 * F1 + iJ12 * F2);
	  dgi = - (iJ21 * F1 + iJ22 * F2);
	  
	  dg = sqrt (dgr*dgr + dgi*dgi);
	  
	  if (dg > smax)
	    {
	      dgr = dgr * smax /dg;
	      dgi = dgi * smax /dg;
	      dg  = smax;
	    }
	  
	  g_r = g_r + lambda * dgr;
	  g_i = g_i + lambda * dgi;
	  
	  SolveLayerEquations ();
	  
	  F1 = real (Deltas) - Delta;
	  F2 = imag (Deltas);
	  
	  F = sqrt (F1*F1 + F2*F2);

	  iter1++;
	}
      while (F > Fold && iter1 < Maxiter);

      lambda = 1.;
      Fold   = F;

      //printf ("g = (%10.3e, %10.3e)  F = (%10.3e, %10.3e) |F| = %10.3e dg = %10.3e\n",
      //      g_r, g_i, F1, F2, F, dg);
    
      iter++;
    }
  while (F > Eta && dg > eps && iter < Maxiter);

  gr = g_r;
  gi = g_i;
}

// ################################
// Target function for zero finding
// ################################
double Layer::Feval (double x)
{
  g_i = x;

  SolveLayerEquations ();

  return imag (Deltas);
}

// ############################################
// Ridder's method for finding root of F(x) = 0
// ############################################
void Layer::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) /2.;    
      double F3 = Feval (x3);
      
      // Iterate x using Ridder's method 
      xold = x;           
      x = x3 - (x3 - x1) * (F2 - F1) * F3 /
	(sqrt (F3 * F3 - F1 * F2) * fabs (F2 - F1));
      Fx = Feval (x);
       
      // Make new value of x upper/lower bound of refined search interval, as appropriate 
      if (Fx * F1 < 0.) 
	{  
	  x2 = x;           
	  F2 = Fx; 
	}
      else 
	{
	  x1 = x;
	  F1 = Fx; 
	}
      iter++;
    } 
  // Iterate until absolute change in x falls below Eta
  while (fabs (x - xold) > Eta && fabs(Fx) > Eta && iter < Maxiter);
}

// ##########################
// Function to read JSON file
// ##########################
json Layer::ReadJSONFile (const string& filename)
{
  ifstream JSONFile (filename);
  json     JSONData;

  if (JSONFile.is_open ())
    {
      try
	{
	  JSONFile >> JSONData;
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

// #####################################
// Function to open new file for writing
// #####################################
FILE* Layer::OpenFilew (const char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL)
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}
