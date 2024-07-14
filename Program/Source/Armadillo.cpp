// Armadillo.cpp

#include "TJ.h"

// Function to find inverse of complex matrix

// ####################################################################################################################
// Function to solve linear system of equations A . X = B, for X, where all quantities are complex rectangular matrices
// ####################################################################################################################
void TJ::SolveLinearSystem (Array<complex<double>,2> A, Array<complex<double>,2> X, Array<complex<double>,2> B)
{
  int size1 = A.extent(0);
  int size2 = A.extent(1);
  int size3 = B.extent(1);

  if (size1 != size2)
    {
      printf ("TJ::SolveLinearSystem: Error - over/underdetermined linear system\n");
      exit (1);
    }
 
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

// ####################################################################################################################
// Function to solve linear system of equations A . X = B, for X, where A is a complex rectangular matrix, and X and B
// are complex vectors
// ####################################################################################################################
void TJ::SolveLinearSystem (Array<complex<double>,2> A, complex<double>* X, complex<double>* B)
{
  int size1 = A.extent(0);
  int size2 = A.extent(1);

  if (size1 != size2)
    {
      printf ("TJ::SolveLinearSystem: Error - over/underdetermined linear system\n");
      exit (1);
    }
 
  // Solve problem using Armadillo
  cx_mat Amat (size1, size1);
  cx_vec Bvec (size1), Xvec (size1);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size1; j++)
      Amat(i, j) = A(i, j);

  for (int i = 0; i < size1; i++)
    Bvec(i) = B[i];
  
  solve (Xvec, Amat, Bvec);

  for (int i = 0; i < size1; i++)
    X[i] = Xvec(i);
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

  cx_mat Bmat = inv (Amat);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      invA(i, j) = Bmat(i, j);
 }

// ###################################################
// Function to return eigenvalues of Hermitian matix H
// ###################################################
void TJ::GetEigenvalues (Array<complex<double>,2> H, double* evals)
{
  int size = H.extent(0);

  // Solve problem using Armadillo
  cx_mat Hmat (size, size);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      Hmat(i, j) = H(i, j);

  vec eigval = eig_sym (Hmat);

  for (int i = 0; i < size; i++)
    evals[i] = eigval(i);
}
