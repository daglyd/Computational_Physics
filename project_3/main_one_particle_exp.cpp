#include "Particle.hpp"
#include "PenningTrap.hpp"


int main(){

    double B0 = 9.65e1;
    double V0 = 9.65e8; 
    double d = 1e4; 

    vec r = vec(3).randn() * 0.1 * d;  // random initial position
    vec v = vec(3).randn() * 0.1 * d;  // random initial velocity
    
    Particle first_particle = Particle(0.99872,40.078 ,r,v); //40.078
         
    PenningTrap trap = PenningTrap(B0,V0,d);
    trap.add_particle(first_particle);

    double time = 100; // ms 
    double n_steps = 10000; 
    double dt = time/n_steps; 

    vec particle_positions = vec(n_steps);
    ofstream file; 
    file.open("Data/main_one_particle_exp/positions_euler_cromer.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;
    for (int i=0; i < n_steps; i++){

        
        file << 0 << ","<< trap.particles.at(0).position_.at(0) 
        <<","<< trap.particles.at(0).position_.at(1)
        <<","<< trap.particles.at(0).position_.at(2) << endl; 
        trap.evolve_forward_Euler(dt);
    }
    file.close();
    // delete &particle_positions;

    file.open("Data/main_one_particle_exp/positions_RK4.csv");
    for (int i=0; i < n_steps; i++){

            file << 0 << "," 
            <<","<< trap.particles.at(0).position_.at(0) 
            <<","<< trap.particles.at(0).position_.at(1)
            <<","<< trap.particles.at(0).position_.at(2) << endl; 
            trap.evolve_RK4(dt);
        }
    

    return 0;
}
