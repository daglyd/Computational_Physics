#include "Particle.hpp"
#include "PenningTrap.hpp"


int main(){
    
    arma_rng::set_seed(12);

    // initial values for the Penniing trap 
    double B0 = 9.65e1;
    double V0 = 9.65e8; 
    double d = 1e4; 
    
    // creating the Penning trap
    PenningTrap trap = PenningTrap(B0,V0,d);
    
    // adding particles to the trap 
    int n_particles = 2; 
    for (int i=0; i<n_particles; i++){
        vec r = vec(3).randn() * 0.01 * d;  // random initial position
        vec v = vec(3).randn() * 0.01 * d;  // random initial velocity

        // Particle particle = Particle(0.99872,40.078 ,r,v); 
        trap.add_particle(Particle(0.99872,40.078 ,r,v));
        cout << "Particle added!" << endl;
    }

    // time and number of steps for the numerical method
    double time = 100; // ms 
    double n_steps = 10000; 
    double dt = time/n_steps; 

    // -------------------
    // Euler-Cromer method
    // ------------------- 
    {
    // opening files for writing the results 
    ofstream file; 
    file.open("Data/main_two_particle_exp/positions_euler_cromer.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;

    ofstream file2; 
    file2.open("Data/main_two_particle_exp/velocities_euler_cromer.csv");
    file2 << "Particle,"<<"vx,"<<"vy,"<<"vz"<<endl;

    // running the numberical method 
    for (int i=0; i < n_steps; i++){

        for (int j=0; j < n_particles; j++){

             // writing positions and velocities to file 
            file << j 
            <<","<<trap.particles.at(j).position_.at(0) 
            <<","<< trap.particles.at(j).position_.at(1)
            <<","<< trap.particles.at(j).position_.at(2) << endl; 
            
            file2 << j 
            <<","<<trap.particles.at(j).velocity_.at(0) 
            <<","<< trap.particles.at(j).velocity_.at(1)
            <<","<< trap.particles.at(j).velocity_.at(2) << endl;       

        }
        // evolves positions and velocites using euler-cromer medthod 
        trap.evolve_Euler_Cromer(dt);
    }
    file.close();
    file2.close();
    }
    // -------------------
    // RungeKutta4 method
    // -------------------
    {
    // opening files for writing the results
    ofstream file; 
    file.open("Data/main_two_particle_exp/positions_RK4.csv");
    file << "Particle,"<<"x,"<<"y,"<<"z"<<endl;

    ofstream file2; 
    file2.open("Data/main_two_particle_exp/velocities_RK4.csv");
    file2 << "Particle,"<<"vx,"<<"vy,"<<"vz"<<endl;

    // running the numberical method 
    for (int i=0; i < n_steps; i++){

        for (int j=0; j < n_particles; j++){

            // writing positions and velocities to file 
            file << j 
            <<","<<trap.particles.at(j).position_.at(0) 
            <<","<< trap.particles.at(j).position_.at(1)
            <<","<< trap.particles.at(j).position_.at(2) << endl; 
            
            file2 << j 
            <<","<<trap.particles.at(j).velocity_.at(0) 
            <<","<< trap.particles.at(j).velocity_.at(1)
            <<","<< trap.particles.at(j).velocity_.at(2) << endl;       

        }
        // evolves positions and velocites using Runge Kutta 4th order method
        trap.evolve_RK4(dt);
    }
    file.close();
    file2.close();
    }
    return 0;
}
