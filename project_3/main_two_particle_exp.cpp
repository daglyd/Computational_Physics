#include "Particle.hpp"
#include "PenningTrap.hpp"


int main(){
    
    arma_rng::set_seed(12);

    double B0 = 9.65e1;
    double V0 = 9.65e8; 
    double d = 1e4; 
    
    PenningTrap trap = PenningTrap(B0,V0,d);
    
    int n_particles = 2; 
    for (int i=0; i<n_particles; i++){
        vec r = vec(3).randn() * 0.01 * d;  // random initial position
        vec v = vec(3).randn() * 0.01 * d;  // random initial velocity

        // Particle particle = Particle(0.99872,40.078 ,r,v); 
        trap.add_particle(Particle(0.99872,40.078 ,r,v));
        cout << "Particle added!" << endl;
    }


    double time = 100; // ms 
    double n_steps = 10000; 
    double dt = time/n_steps; 

    // Euler-Cromer method 
    {
    ofstream file; 
    file.open("Data/main_two_particle_exp/positions_euler_cromer.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;

    ofstream file2; 
    file2.open("Data/main_two_particle_exp/velocities_euler_cromer.csv");
    file2 << "Particle,"<<"vx,"<<"vy,"<<"vz"<<endl;

    for (int i=0; i < n_steps; i++){

        for (int j=0; j < n_particles; j++){
            file << j 
            <<","<<trap.particles.at(j).position_.at(0) 
            <<","<< trap.particles.at(j).position_.at(1)
            <<","<< trap.particles.at(j).position_.at(2) << endl; 
            
            file2 << j 
            <<","<<trap.particles.at(j).velocity_.at(0) 
            <<","<< trap.particles.at(j).velocity_.at(1)
            <<","<< trap.particles.at(j).velocity_.at(2) << endl;       

        }
        
        trap.evolve_Euler_Cromer(dt);
    }
    file.close();
    file2.close();
    }
    // RungeKutta4 method
    {
    ofstream file; 
    file.open("Data/main_two_particle_exp/positions_RK4.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;

    ofstream file2; 
    file2.open("Data/main_two_particle_exp/velocities_RK4.csv");
    file2 << "Particle,"<<"vx,"<<"vy,"<<"vz"<<endl;

    for (int i=0; i < n_steps; i++){

        for (int j=0; j < n_particles; j++){
            file << j 
            <<","<<trap.particles.at(j).position_.at(0) 
            <<","<< trap.particles.at(j).position_.at(1)
            <<","<< trap.particles.at(j).position_.at(2) << endl; 
            
            file2 << j 
            <<","<<trap.particles.at(j).velocity_.at(0) 
            <<","<< trap.particles.at(j).velocity_.at(1)
            <<","<< trap.particles.at(j).velocity_.at(2) << endl;       

        }
        
        trap.evolve_RK4(dt);
    }
    file.close();
    file2.close();
    }
    return 0;
}
