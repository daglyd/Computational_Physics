#include "PenningTrap.hpp"


/// <summary>Penning trap object for modelling a Penning trap</summary>
/// <param name="B0_in">double: Magnitude of the magnetic field in the trap</param>  
/// <param name="V0_in">double: Magnitude of the electric field</param>  
/// <param name="d_in">double: Characteristic dimension of the trap</param> 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{   
    
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
}

/// <summary> adds a particle to the trap </summary>
/// <param name="p_in">Particle: particle to be added</param> 
void PenningTrap::add_particle(Particle p_in){
    particles.push_back(p_in);
}

/// <summary> Counts the number of particles in the trap </summary>
/// <returns>int: number of particles in trap</returns>  
int PenningTrap::particles_count()
{
    return particles.size();
}

/// <summary> Returns the r-vector </summary>
vec PenningTrap::return_r(){
    vec r = vec(particles.size());
    int n = r.size();
    for (int i=0;i<n;i++){
    }
    return r;
}

/// <summary> Evaluates the external E-field at r(x,y,z) </summary>
/// <param name="r">arma::vec: vector of positions r</param> 
/// <returns>arma::vec: E-field as vector </returns>  
vec PenningTrap::external_E_field(vec r)
{
    
    vec c = {-2,-2,4};
    // vec E = vec(3, fill::zeros);

    vec E = (c%r) * ( V0_ / (2*pow(d_,2)) );
    
    return E;
}

/// <summary> Evaluates the external B-field at r(x,y,z) </summary>
/// <param name="r">arma::vec: vector of positions r</param> 
/// <returns>arma::vec: B-field as vector </returns>  
vec PenningTrap::external_B_field(vec r){

    vec B = vec(r.size(), fill::zeros);


    B.at(3-1) = B0_;
    return B; 
}


/// <summary> Evaluates the force on particle i from j </summary>
/// <param name="i">int: particle nr i </param> 
/// <param name="j">int: particle nr j </param> 
/// <returns>arma::vec: Force from particle j on i as vector </returns> 
vec PenningTrap::force_particle(int i, int j){

    Particle particle_i = particles.at(i);
    Particle particle_j = particles.at(j);

    vec E = vec(3);

    E  = k_e *  (particle_j.charge_ * 
        (particle_i.position_ - particle_j.position_) /
        pow(norm(particle_i.position_ - particle_j.position_),3));

    return particle_i.charge_ * E; 
}

/// <summary> Evaluates the total external forces on particle i </summary>
/// <param name="i">int: particle nr i </param> 
/// <returns>arma::vec: Total external force on particle i as vecotr </returns> 
vec PenningTrap::total_force_external(int i){

    Particle particle_i = particles.at(i);

    vec a = particle_i.charge_*external_E_field(particle_i.position_);
    vec b = particle_i.charge_*particle_i.velocity_;
    vec c = external_B_field(particle_i.position_);
    vec force = a + cross(b,c);
    return force; 
}


/// <summary> Evaluates the total forces on particle i from other particles </summary>
/// <param name="i">int: particle nr i </param> 
/// <returns>arma::vec: Total force on particle i from other particles as vector </returns>
vec PenningTrap::total_force_particles(int i){
    Particle particle_i = particles.at(i);
    vec E = vec(3,fill::zeros);

    for (int j=0; j < particles.size(); j++){
        if (j == i){
            ;
        }
        else{
        E += particles.at(j).charge_ * 
            (particle_i.position_ - particles.at(j).position_) /
            pow(norm(particle_i.position_ - particles.at(j).position_),3);
        }
    }
    
    return particle_i.charge_ * E;
}

/// <summary> Evaluates the total forces on particle i from other particles </summary>
/// <param name="i">int: particle nr i </param> 
/// <returns>arma::vec: Total force on particle i from other particles as vector </returns>
vec PenningTrap::total_force(int i){

    vec force = total_force_external(i) + total_force_particles(i);
    return force; 
}

/// <summary> Evolves the particle positions and velocities one step using the forward-Euler method </summary>
/// <param name="dt">double: time step dt </param> 
void PenningTrap::evolve_forward_Euler(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){

        vec a = total_force(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}


/// <summary> Evolves the particle positions and velocities one step using the Euler-Cromer method </summary>
/// <param name="dt">double: time step dt </param> 
void PenningTrap::evolve_Euler_Cromer(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){

        vec a = total_force(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}


/// <summary> Evolves the particle positions and velocities one step using the Euler-Cromer method 
/// using only external forces </summary>
/// <param name="dt">double: time step dt </param> 
void PenningTrap::evolve_Euler_Cromer_without_particle_interaction(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){
    //    Particle p = particles.at(i);

        vec a = total_force_external(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}


/// <summary> Evolves the particle positions and velocities one step using the Runge Kutta 4th order method </summary>
/// <param name="dt">double: time step dt </param> 
void PenningTrap::evolve_RK4(double dt){
    
    int n = particles.size();
    double half_dt = dt*0.5;
    for (int i=0; i< n; i++){
        
        vec velocity_n = particles.at(i).velocity_;
        vec position_n = particles.at(i).position_; 

        vec k1 = total_force(i)/particles.at(i).mass_;

        particles.at(i).velocity_ = velocity_n + k1 * half_dt;
        particles.at(i).position_ = position_n + k1 * half_dt;
        vec k2 = total_force(i)/particles.at(i).mass_;
        
        particles.at(i).velocity_ = velocity_n + k2 * half_dt;
        particles.at(i).position_ = position_n + k2 * half_dt;     
        vec k3 = total_force(i)/particles.at(i).mass_;

        particles.at(i).velocity_ = velocity_n + k3 * dt;
        particles.at(i).position_ = position_n + k3 * dt;    
        vec k4 = total_force(i)/particles.at(i).mass_;

        particles.at(i).velocity_ = velocity_n + (dt/6) * (k1 + 2*(k2+k3) + k4); 
        particles.at(i).position_ = position_n + (dt/6) * (k1 + 2*(k2+k3) + k4); 

    }
}