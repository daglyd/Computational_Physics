#include <iostream>
#include "thomas_algo.hpp"
#include <string>
#include <chrono>

using namespace std;
using namespace arma;



int main(){
    
    // init siganture
    int signature[3] = {-1,2,-1};
    
    // vectors for recoding time and n steps 
    vec times_1 = vec(6*10);
    vec steps_1 = vec(6*10);
    int j = 0; 
    for (int s=0; s<10;s++){
        
        // running general algorithm 10 times
        // recoding times and steps 
        for (int i=10; i <= 1e6; i *= 10){
            times_1[j] = thomas_algo_general(i, signature);
            steps_1[j] = i;
            j ++;
        }
    }

    // saving steps and times to csv file 
    mat steps_times_1 = join_rows(steps_1,times_1);
    string filename_1 = "Data/general_algo_times.csv";
    ofstream file(filename_1);
    file << "steps,times" << endl;
    steps_times_1.save(file,csv_ascii);
    

    // vectors for recording time and n steps for special algorithm
    vec times_2 = vec(6*10);
    vec steps_2 = vec(6*10);

    int k = 0; 
    for (int s=0; s<10; s++){
    for (int i=10; i <= 1e6; i *= 10){

        // running special algorithm 10 times
        // recoding times and steps 
        times_2[k] = thomas_algo_special(i, signature);
        steps_2[k] = i;
        k ++;
    }
    }

    // saving steps and times to csv file 
    mat steps_times_2 = join_rows(steps_2,times_2);
    string filename_2 = "Data/special_algo_times.csv";
    ofstream file_2(filename_2);
    file_2 << "steps,times" << endl;
    steps_times_2.save(file_2,csv_ascii);

    return 0;
}


