
#ifndef __thomas_algo_hpp__
#define __thomas_algo_hpp__

#include <iostream>
#include <math.h>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;

double thomas_algo_general(int n, int *signature);
double thomas_algo_special(int n, int *signature);

mat forward_sub_special(vec a, vec b, vec c, vec g);
mat backward_sub_special(vec b_, vec g_, vec c);
mat forward_sub(vec a, vec b, vec c, vec g);
mat backward_sub(vec b, vec g, vec c);

mat tridiagonal_matrix(int *signature, int n);
vec f(vec x, double h);

#endif