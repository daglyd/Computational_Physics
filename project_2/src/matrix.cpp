#include "headers.hpp"

double max_offdiag_symmetric(const mat& A, int& k, int &l){

    auto t1 = chrono::high_resolution_clock::now();
    
    mat U_matrix = trimatu(A,1);
    U_matrix = abs(U_matrix);
    

    double max_value = abs(U_matrix).max();
    uvec max_value_index = find(U_matrix == max_value,k=1);

    int index = conv_to<int>::from(max_value_index);

    k = index % U_matrix.n_rows;
    l = index / U_matrix.n_rows;


    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

    return U_matrix[index];
}

mat tridiagonal_matrix(double *signature, int n){

    // init diagonal vectors of matrix A, (a,b,c)
    // vec a = vec(n-1);
    // vec b = vec(n);
    // vec c = vec(n-1);

    // // filling vectors according to signature
    // for (int i=0; i<n; i++){
    //     a[i+1] = signature[0];
    //     b[i] = signature[1];
    //     c[i+1] = signature[2];
    // }
    
    mat matrix = mat(n,n, fill::zeros);
    matrix.diag() += signature[1];
    matrix.diag(1) += signature[0];
    matrix.diag(-1) += signature[2];


    // setting first element of a to 0 and last elem of c.
    // this to equate lenght of a & c with b. 
    // a.head(1) = 0;
    // c.tail(1) = 0;
    
    // // returning joined matrix 
    // mat m = join_rows(a,b);
    // m = join_rows(m,c);
    return matrix;

}