// Integrate.cpp

// PROGRAM ORGANIZATION:
//
// void Flux:: CalcQGP                  ()
// void Flux:: CalcrP                   ()
// void Flux:: CalcStraightAngleGeneral ()
// int  Flux:: Rhs1                     (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs2                     (double r, const double y[], double dydr[], void*)
// int  Flux:: Rhs3                     (double r, const double y[], double dydr[], void*)

#include "Flux.h"

// #######################################
// Function to calculate q(P)/g(P) profile
// #######################################
void Flux::CalcQGP ()
{
  gsl_odeiv2_system           sys1 = {pRhs1, NULL, 4, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys1, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [4]; 
  double                      r;

  for (int j = 0; j < NPSI; j++)
    {
      r    = 0.;
      y[0] = RP[j];
      y[1] = Zaxis;
      y[2] = 0.;
      y[3] = 0.;
      qgp  = GP[j];
  
      int status = gsl_odeiv2_driver_apply (d, &r, 2.*M_PI, y);

      if (status != GSL_SUCCESS)
	{
	  printf ("CalcGGP: status != GSL_SUCCESS\n");
	  exit (1);
	}
      
      QGP[j] = y[2];
      QP [j] = QGP[j]*GP[j];

      if (j%10 == 0 || j > NPSI-10 || j < 10)
	{
	  printf ("j = %4d  PsiN = %11.4e  q = %11.4e\n", j, 1.-P[j], QP[j]);
	  fflush (stdout);
	}
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ##################################
// Function to calculate r(P) profile
// ##################################
void Flux::CalcrP ()
{
  gsl_odeiv2_system           sys2 = {pRhs2, NULL, 1, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys2, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [1]; 
  double                      r;

  // Calculate r[Psi]
  for (int j = 0; j < NPSI; j++)
    {
      r    = P[j];
      y[0] = 0.;
      
      int status = gsl_odeiv2_driver_apply (d, &r, 1., y);
	  
      if (status != GSL_SUCCESS)
	{
	  printf ("CalcrP: status != GSL_SUCCESS\n");
	  exit (1);
	}
      
      rP[j] = sqrt (y[0]);
    }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// ###########################################################################
// Function to calculate straight angle data at general magnetic-flux surfaces
// ###########################################################################
void Flux::CalcStraightAngleGeneral ()
{
  gsl_odeiv2_system           sys3 = {pRhs3, NULL, 2, this};
  const gsl_odeiv2_step_type* T    = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver*          d    = gsl_odeiv2_driver_alloc_y_new (&sys3, T, H0, ACC/2., ACC/2.);
  double*                     y    = new double [2]; 
  double                      r;
  
  for (int j = 0; j < NPSI; j++)
    {
      r    = 0.;
      y[0] = RP[j];
      y[1] = Zaxis;
      qgp  = fabs (Psic) * QP[j] /GP[j];
	      
      gsl_matrix_set (RRst, j, 0, y[0]);
      gsl_matrix_set (ZZst, j, 0, y[1]);

      for (int k = 1; k < NTHETA; k++)
	{
	  int status = gsl_odeiv2_driver_apply (d, &r, th[k], y); 
	  
	  if (status != GSL_SUCCESS)
	    {
	      printf ("CalcStraightAngleGeneral: status != GSL_SUCCESS\n");
		  exit (1);
	    }

	  gsl_matrix_set (RRst, j, k, y[0]);
	  gsl_matrix_set (ZZst, j, k, y[1]);
	}

      if (j%10 == 0 || j > NPSI-10 | j < 10)
	{
	  printf ("j = %4d  PsiN = %11.4e  r/a = %11.4e\n", j, 1.-P[j], rP[j] /ra);
	  fflush (stdout);
	}
     }

  gsl_odeiv2_driver_free (d);
  delete[]                y;
}

// #####################################################
// Function to evaluate right-hand sides of q/g equation
// #####################################################
int Flux::Rhs1 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - q/g
  // y[3] - 1/gamma
  
  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Grad = PsiR*PsiR + PsiZ*PsiZ;
  double RB   = sqrt (qgp*qgp + Psic*Psic * Grad);
  double Rc   = Raxis - y[0];
  double Zc   = y[1]  - Zaxis;
  double fac  = (Rc*Rc + Zc*Zc) /(- Zc*PsiZ + Rc*PsiR);

  dydr[0] = - PsiZ * fac;
  dydr[1] = + PsiR * fac;
  dydr[2] = (1./2./M_PI /fabs(Psic)) * fac /y[0];
  dydr[3] = (1./2./M_PI /fabs(Psic)) * RB * fac;
    
  return GSL_SUCCESS;
}

// ###################################################
// Function to evaluate right-hand sides of r equation
// ###################################################
int Flux::Rhs2 (double r, const double y[], double dydr[], void*)
{
  // y[0] - (rP)^2    

  double qg = Interpolate (NPSI, PsiN, QGP, 1.-r, 0);

  dydr[0] = 2.*fabs(Psic) * qg;

  return GSL_SUCCESS;
}

// #######################################################
// Function to evaluate right-hand sides of theta equation
// #######################################################
int Flux::Rhs3 (double r, const double y[], double dydr[], void*)
{
  // y[0] - R     
  // y[1] - Z

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);

  dydr[0] = - qgp * y[0] * PsiZ;
  dydr[1] = + qgp * y[0] * PsiR;

  return GSL_SUCCESS;
}

