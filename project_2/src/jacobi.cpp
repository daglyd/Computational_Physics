#include "headers.hpp"


void jacobi_eigensolver( mat& A, double eps, vec& eigenvalues, mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    mat R = mat(A.n_rows,A.n_rows, fill::eye);

    int k,l;
    

    double max_value = max_offdiag_symmetric(A,k,l);

    while ((max_value > eps) && (iterations <= maxiter)){
        jacobi_rotate(A, R, k, l);
        max_value = max_offdiag_symmetric(A,k,l);
        
        iterations ++;
        }
    if (max_value < eps){
        converged = true;
    }
    else{
        converged = false;
    }
   

    eigenvalues = A.diag();
    eigenvectors = R; 


}


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int  l){
    double t,c,s;

    if (A.at(k,l)!=0.0){
        double tau = (A.at(l,l) - A.at(k,k))/(2.0*A.at(k,l));

        if (tau > 0.0){
            t = 1.0/(tau + sqrt(1.0 + pow(tau,2))); // -tau + sqrt(1+pow(tau,2));
        }
        else{
            t = -1.0/( -tau + sqrt(1.0 + pow(tau,2))); //-tau - sqrt(1+pow(tau,2));
        }
        c = 1.0 / sqrt(1.0+pow(t,2.0));
        s = c*t;
        }
        else{
            c = 1.0;
            s = 0.0;
            t = 0.0;
        }


        A.at(k,k) = A.at(k,k)*pow(c,2.0) - 2.0*A.at(k,l)*c*s + A.at(l,l)*pow(s,2.0);
        A.at(l,l) = A.at(l,l)*pow(c,2.0) + 2.0*A.at(k,l)*c*s + A.at(k,k)*pow(s,2.0);

        A.at(k,l) = 0.0;
        A.at(l,k) = 0.0;

        for (int i = 0; i< A.n_rows; i++){
            if (i!=k && i!=l){
                double A_ik = A.at(i,k);
                double A_il = A.at(i,l);

                A.at(i,k) = A_ik*c - A_il*s;
                A.at(k,i) = A.at(i,k);

                A.at(i,l) = A_il*c + A_ik*s;
                A.at(l,i) = A.at(i,l);
            }
            double R_ik = R.at(i,k);
            double R_il = R.at(i,l);
            R.at(i,k) = R_ik*c - R_il*s;
            R.at(i,l) = R_il*c + R_ik*s;

        }
}