#include "Particle.hpp"
#include "PenningTrap.hpp"


int main(){

    // initial values for the Penniing trap 
    double B0 = 9.65e1;
    double V0 = 9.65e8; 
    double d = 1e4; 

    // random generated initial positions and velocities 
    vec r = vec(3).randn() * 0.1 * d;  // random initial position
    vec v = vec(3).randn() * 0.1 * d;  // random initial velocity
    
    // first particle object 
    Particle first_particle = Particle(0.99872,40.078 ,r,v); //40.078
         
    // creating the Penning trap
    PenningTrap trap = PenningTrap(B0,V0,d);

    // add particle to trap 
    trap.add_particle(first_particle);

    // time and number of steps for the numerical method
    double time = 100; // ms 
    double n_steps = 10000; 
    double dt = time/n_steps; 

    // opening file for writing the results 
    vec particle_positions = vec(n_steps);
    ofstream file; 
    file.open("Data/main_one_particle_exp/positions_euler_cromer.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;

    // running the numberical method 
    for (int i=0; i < n_steps; i++){

        // writing positions and velocities to file 
        file << 0 << ","<< trap.particles.at(0).position_.at(0) 
        <<","<< trap.particles.at(0).position_.at(1)
        <<","<< trap.particles.at(0).position_.at(2) << endl; 
        
        // evolves positions and velocites using forward euler
        trap.evolve_forward_Euler(dt);
    }
    file.close();

    // opening file for writing the results 
    file.open("Data/main_one_particle_exp/positions_RK4.csv");
    for (int i=0; i < n_steps; i++){

            // writing positions and velocities to file 
            file << 0 << "," 
            <<","<< trap.particles.at(0).position_.at(0) 
            <<","<< trap.particles.at(0).position_.at(1)
            <<","<< trap.particles.at(0).position_.at(2) << endl; 

            // evolves positions and velocites using Runge Kutta 4th order method
            trap.evolve_RK4(dt);
        }
    

    return 0;
}
