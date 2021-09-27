#include "src/headers.hpp"

int main(){

    cout << "Problem 3" << endl;
    cout << "-----------------"<< endl;
    problem_3();

    cout << "===============================================" << endl;

    cout << "Problem 4" << endl;
    cout << "-----------------"<< endl;
    problem_4();

    cout << "===============================================" << endl;

    cout << "Problem 5" << endl;
    cout << "-----------------"<< endl;
    problem_5();

    cout << "===============================================" << endl;

    cout << "Problem 6" << endl;
    cout << "-----------------"<< endl;

    int n = 15;
    mat iterations = mat(n,2);
    int j = 5;
    for (int i=0;i<=n; i++){
        iterations.at(i,1) = problem_6(j);
        iterations.at(i,0) = j;
        j +=5;
    }
    string filename = "Data/problem_6/iterations_N.csv";
    ofstream file(filename);
    file << "N,iterations" << endl;
    iterations.save(file, csv_ascii);

    cout << "===============================================" << endl;

    cout << "Problem 7" << endl;
    cout << "-----------------"<< endl;
    problem_7(10);
    problem_7(100);
    cout << "Saved to file in directory ./Data/" << endl;

    return 0;
}