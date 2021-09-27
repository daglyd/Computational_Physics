#ifndef __headers_hpp__
#define __headers_hpp__

#include <iostream>
#include <math.h>
#include <assert.h>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;

mat tridiagonal_matrix(double *signature, int n);

int analytical(vec& eigval, mat& eigvec);
int problem_6(double n);

double max_offdiag_symmetric(const mat&, int& k, int &l);


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver( mat& A, double eps, vec& eigenvalues, mat& eigenvectors, const int maxiter, int& iterations, bool& converged);
void problem_3();
void problem_5();
void problem_4();
void problem_7(int x);


#endif