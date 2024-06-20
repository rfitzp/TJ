// Armadillo.cpp

#include "TJ.h"

// Function to find inverse of complex matrix

// ####################################################################################################################
// Function to solve linear system of equations A . X = B, for B, where all quantities are complex rectangular matrices
// ####################################################################################################################
void TJ::SolveLinearSystem (Array<complex<double>,2> A, Array<complex<double>,2> X, Array<complex<double>,2> B)
{
  int size1 = A.extent(0);
  int size2 = A.extent(1);
  int size3 = B.extent(1);
 
  // Solve problem using Armadillo
  cx_mat Amat (size1, size2), Xmat (size2, size3), Bmat (size1, size3);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      Amat(i, j) = A(i, j);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size3; j++)
      Bmat(i, j) = B(i, j);
  
  solve (Xmat, Amat, Bmat);

  for (int i = 0; i < size2; i++)
    for (int j = 0; j < size3; j++)
      X(i, j) = Xmat(i, j);
 }

// ################################
// Function to invert square matrix
// ################################
void TJ::InvertMatrix (Array<complex<double>,2> A, Array<complex<double>,2> invA)
{
  int size = A.extent(0);
 
  // Solve problem using Armadillo
  cx_mat Amat (size, size);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      Amat(i, j) = A(i, j);

  cx_mat Bmat = inv(Amat);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      invA(i, j) = Bmat(i, j);
 }

