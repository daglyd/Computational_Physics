#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

double u(double x);


int main()
  {
    int n = 1000;
    vec x = linspace(0,1,n);
    vec y = vec(n);

    for (int i=0; i<x.n_elem; i++){
      y[i] = u(x[i]);
    }

    mat z = join_rows(x,y);
    
    ofstream file("Data/prob_2_data.csv");
    file << "x,u(x)" << endl;
    z.save(file, csv_ascii);

    return 0;
  }

double u(double x){
    double u = 1 - (1 - exp(-10))*x - exp(-10*x);
    return u; 
}