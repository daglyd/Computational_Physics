#include "thomas_algo.hpp"

double thomas_algo_special(int n, int *signature){

    // creating vector x with n elements and calculating step size
    vec x = linspace(0,1, n);
    double h = (x.tail(1)[0] - x.head(1)[0]) / n;

    // initializing g(x) and tridiagonal matrix A
    vec g = f(x,h);
    mat m = tridiagonal_matrix(signature, n);

    // Special: running the forward and backward substituion & timing the algorithm
    auto t1 = chrono::high_resolution_clock::now();
    mat b_g_mat = forward_sub_special(m.col(0), m.col(1), m.col(2), g);
    vec v =  backward_sub_special(b_g_mat.col(0), b_g_mat.col(1), m.col(2));
    auto t2 = chrono::high_resolution_clock::now();

    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

    // joining the results v(x) with x to a matrix and writing to csv file
    mat z = join_rows(x,v);
    string filename = "Data/prob9_"+to_string(n)+"_data.csv";
    ofstream file(filename);
    file << "x,v(x)" << endl;
    z.save(file, csv_ascii);

    // returning the duration of the algorithm 
    return duration_seconds;
}

double thomas_algo_general(int n, int *signature ){
    
    // creating vector x with n elements and calculating step size
    vec x = linspace(0,1, n);
    double h = (x.tail(1)[0] - x.head(1)[0]) / n;
    
    // initializing g(x) and tridiagonal matrix A
    vec g = f(x,h);
    mat m = tridiagonal_matrix(signature, n);

    // General: running the forward and backward substituion & timing the algorithm
    auto t1 = chrono::high_resolution_clock::now();
    mat b_g_mat = forward_sub(m.col(0), m.col(1), m.col(2), g);
    vec v =  backward_sub(b_g_mat.col(0), b_g_mat.col(1), m.col(2));
    auto t2 = chrono::high_resolution_clock::now();

    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

    // joining the results v(x) with x to a matrix and writing to csv file
    mat z = join_rows(x,v);
    string filename = "Data/prob7_"+to_string(n)+"_data.csv";
    ofstream file(filename);
    file << "x,v(x)" << endl;
    z.save(file, csv_ascii);

    // returning the duration of the algorithm
    return duration_seconds;
}

mat forward_sub_special(vec a, vec b, vec c, vec g){

    // init. the new vectors for b and g 
    vec b_ = vec(b.n_elem);
    vec g_ = vec(g.n_elem);
    b_[0] = b[0];
    g_[0] = g[0];

    // running the special forward subs. 
    for (int i = 1; i <b.n_elem; i++){
        b_[i] = b[i] - (1/b_[i-1]);
        g_[i] = g[i] + (g_[i-1]/b_[i-1]);
    }

    // joining new b and g vectors to single matrix and returning matrix
    mat m = join_rows(b_,g_);
    return m;
}
mat forward_sub(vec a, vec b, vec c, vec g){

    // init. the new vectors for b and g 
    vec b_ = vec(b.n_elem);
    vec g_ = vec(g.n_elem);
    b_[0] = b[0];
    g_[0] = g[0];

    // running the general forward subs.
    for (int i = 1; i <b.n_elem; i++){
        b_[i] = b[i] - (a[i]/b_[i-1])*c[i-1];
        g_[i] = g[i] - (a[i]/b_[i-1])*g_[i-1];
    }

    // joining new b and g vectors to single matrix and returning matrix
    mat m = join_rows(b_,g_);
    return m;
}

mat backward_sub_special(vec b_, vec g_, vec c){

    // init v vector 
    vec v = vec(g_.n_elem);

    // init last item in v 
    v.tail(1) = g_.tail(1)/b_.tail(1);

    // running special backward subs. 
    for (int i = g_.n_elem -1; i >=0; i--){
        v[i] = (g_[i] + v[i+1])/b_[i];
    }
    
    return v;
}

mat backward_sub(vec b_, vec g_, vec c){

    // init v vector
    vec v = vec(g_.n_elem);

    // init last item in v 
    v.tail(1) = g_.tail(1)/b_.tail(1);

    // running general backward subs. 
    for (int i = g_.n_elem -1; i >=0; i--){
        v[i] = (g_[i] - c[i]*v[i+1])/b_[i];
    }
    
    
    return v;
}

vec f(vec x, double h){
    
    // init righ hand side function 
    vec f = vec(x.n_elem);

    // filling f vector -> 100*e^{-10x} * h^2
    for (int i=0; i<x.n_elem; i++){
        f[i] = 100*exp(-10*x[i])*pow(h,2);
    }

    return f;
}


mat tridiagonal_matrix(int *signature, int n){

    // init diagonal vectors of matrix A, (a,b,c)
    vec a = vec(n);
    vec b = vec(n);
    vec c = vec(n);

    // filling vectors according to signature
    for (int i=0; i<n; i++){
        a[i] = signature[0];
        b[i] = signature[1];
        c[i] = signature[2];
    }

    // setting first element of a to 0 and last elem of c.
    // this to equate lenght of a & c with b. 
    a.head(1) = 0;
    c.tail(1) = 0;
    
    // returning joined matrix 
    mat m = join_rows(a,b);
    m = join_rows(m,c);
    return m;

}