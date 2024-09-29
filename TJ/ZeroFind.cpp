// ZeroFind.cpp

#include "TJ.h"

// ####################################################################
//  Routine to find approximate root of F(x) = 0 using Ridder's method. 
//  Search takes place in interval (0., 1.). 
//  Interval is chopped into nint equal segments.
// ####################################################################
double TJ::RootFind ()
{
  double x1 = EPS, x2 = 1., F1, F2 = 0., root;

  // Chop search interval into Nr segments  
  for (int i = 0; i < Nr; i++)
    {
      double x1_seg = x1 + (x2 - x1) * double (i)     / double (Nr);
      double x2_seg = x1 + (x2 - x1) * double (i + 1) / double (Nr);
      
      if (i == 0) 
	F1 = Feval (x1_seg);
      else 
	F1 = F2;
      F2 = Feval (x2_seg);
      
      // Call Ridder's method for segment containing zero
      if (F1 * F2 < 0.)
	Ridder (x1_seg, x2_seg, F1, F2, root);
    }

  return root;
}

// #############################################
// Ridder's method for finding root of F(x) = 0.
// #############################################
void TJ::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) / 2.;    
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


