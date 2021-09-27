#include "headers.hpp"



void problem_7(int x){
    double n = x -1;
    double h = 1.0/n; 
    double signature[3] = {-1.0/pow(h,2.0),2.0/pow(h,2.0),-1.0/pow(h,2.0)}; 

    mat A = tridiagonal_matrix(signature,n);

    double eps = 1e-8;
    vec eigenvalues = vec(n);
    mat eigenvectors = mat(n,n);

    int maxiter = 10000;
    int iterations; 
    bool converged;


    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);


    vec norm_eigval = normalise(eigenvalues);
    mat norm_eigvec = normalise(eigenvectors);

    string filename = "Data/problem_7/Eigenvalues_"+to_string(int(n))+".csv";
    ofstream file(filename);
    norm_eigval.save(file, csv_ascii);

    string filename2 = "Data/problem_7/Eigenvectors_"+to_string(int(n))+".csv";
    ofstream file2(filename2);
    norm_eigvec.save(file2, csv_ascii); 

    
}

int problem_6(double n){
    // double n = n;
    double h = 1.0/n; 
    double signature[3] = {-1.0/pow(h,2.0),2.0/pow(h,2.0),-1.0/pow(h,2.0)}; 

    mat A = tridiagonal_matrix(signature,n);

    double eps = 1e-8;
    vec eigenvalues = vec(n);
    mat eigenvectors = mat(n,n);

    int maxiter = 100000;
    int iterations; 
    bool converged;


    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    return iterations;
}
void problem_5(){

    double n = 6.0;
    double h = 1.0/n; 
    double signature[3] = {-1.0/pow(h,2.0),2.0/pow(h,2.0),-1.0/pow(h,2.0)}; 

    mat A = tridiagonal_matrix(signature,n);
    double eps = 1e-8;
    vec eigenvalues = vec(n);
    mat eigenvectors = mat(n,n);
    vec eigval_analyt(n);
    mat eigvec_analyt(n,n);

    int maxiter = 100;
    int iterations; 
    bool converged;

    cout << "Matrix A ("+to_string(int(n))+"X"+to_string(int(n))+")" << endl;
    cout << "-----------------"<<  endl;
    cout << A << endl;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    analytical(eigval_analyt, eigvec_analyt);
    cout << "Eigenvalues - Analytical Eigenvalues" << endl;
    cout << "------------------------------------"<< endl;
    cout << join_rows(normalise(sort(eigenvalues)), normalise(eigval_analyt)) << endl;

    cout << "Eigenvectors" << endl;
    cout << "------------"<<  endl;
    cout << normalise(sort(eigenvectors)) << endl;

    cout << "Analytical Eigenvectors" << endl;
    cout << "-----------------------"<< endl;
    cout << normalise(sort(eigvec_analyt)) << endl;

    string filename = "Data/problem_5/Eigenvalues.csv";
    ofstream file(filename);
    eigenvalues.save(file, csv_ascii);

    string filename2 = "Data/problem_5/Eigenvectors.csv";
    ofstream file2(filename2);
    eigenvectors.save(file2, csv_ascii);
}
void problem_4(){

    int k,l;

     mat matrix = {1,0,0,0.5,
                0,1,-0.7,0,
                0,-0.7,1,0,
                0.5,0,0,1};
    matrix.reshape(4,4);
    matrix = trans(matrix);

    double max_value = max_offdiag_symmetric(matrix,k,l);

    assert(k==1);
    assert(l==2);
    assert(max_value==0.7);
    cout << "Test of max_offdiag_symmetric() ran successfully" << endl;
}


void problem_3(){
    double n = 6;
    double h = 1/n; 
    double signature[3] = {-1/pow(h,2),2/pow(h,2),-1/pow(h,2)}; 

    mat matrix = tridiagonal_matrix(signature, n);
    vec eigval;
    mat eigvec;
    vec eigval_analyt(n);
    mat eigvec_analyt(n,n);
    
    cout << "6x6 Matrix:" << endl;
    cout << "------------"<< endl;
    cout << matrix << endl;

    eig_sym(eigval, eigvec, matrix);

    analytical(eigval_analyt, eigvec_analyt);
    cout << "Eigenvalues - Analytical Eigenvalues" << endl;
    cout << "------------------------------------"<< endl;
    cout << join_rows(normalise(eigval), normalise(eigval_analyt)) << endl;


    cout << "Eigenvectors" << endl;
    cout << "------------"<< endl;
    cout << normalise(eigvec) << endl;

    cout << "Analytical Eigenvectors" << endl;
    cout << "-----------------------"<< endl;
    cout << normalise(eigvec_analyt) << endl;

   

}

int analytical(vec& eigval, mat& eigvec){
    
    double n = 6.0;
    double h = 1/n;
    double a = -1.0/pow(h,2.0);
    double d = 2.0/pow(h,2.0);
    
    for (int i = 1; i<=n; i++){
        eigval.at(i-1) = d + 2.0*a*cos( (i*M_PI)/(n+1) );

        for (int j = 1; j<=n; j++){
            eigvec.at(i-1,j-1) = sin((j*i*M_PI) / (n+1) );    
        }
    }
    return 0;
}

