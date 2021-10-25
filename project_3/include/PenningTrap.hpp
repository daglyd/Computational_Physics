#ifndef __PenningTrap__hpp__
#define __PenningTrap__hpp__

#include <armadillo>
#include <vector>
#include <math.h>
#include "Particle.hpp"

using namespace std;
using namespace arma; 

class PenningTrap
{
    private:
     // No private variables
    // vector<Particle> particle_collection;

    public:

        double B0_, V0_, d_; 
        vector<Particle> particles;
        double k_e = 1.38935333e5;

        // Constructor
        PenningTrap(double B0_in, double V0_in, double d_in);

        // Add a particle to the trap
        void add_particle(Particle p_in);

        int particles_count();

        vec return_r();

        // External electric field at point r=(x,y,z)
        arma::vec external_E_field(arma::vec r);  

        // External magnetic field at point r=(x,y,z)
        arma::vec external_B_field(arma::vec r);  

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields
        arma::vec total_force_external(int i);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec total_force(int i);

        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt);

        void evolve_Euler_Cromer(double dt);

        void evolve_Euler_Cromer_without_particle_interaction(double dt);
};

#endif