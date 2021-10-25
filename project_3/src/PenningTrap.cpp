#include "PenningTrap.hpp"


// Penning trap class
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{   
    
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
}

void PenningTrap::add_particle(Particle p_in){
    particles.push_back(p_in);
}

int PenningTrap::particles_count()
{
    return particles.size();
}

vec PenningTrap::return_r(){
    vec r = vec(particles.size());
    int n = r.size();
    for (int i=0;i<n;i++){
    }
    return r;
}
vec PenningTrap::external_E_field(vec r)
{
    
    vec c = {-2,-2,4};
    // vec E = vec(3, fill::zeros);

    vec E = (c%r) * ( V0_ / (2*pow(d_,2)) );
    
    return E;
}

vec PenningTrap::external_B_field(vec r){

    vec B = vec(r.size(), fill::zeros);


    B.at(3-1) = B0_;
    return B; 
}

vec PenningTrap::force_particle(int i, int j){

    Particle particle_i = particles.at(i);
    Particle particle_j = particles.at(j);

    vec E = vec(3);

    E  = k_e *  (particle_j.charge_ * 
        (particle_i.position_ - particle_j.position_) /
        pow(norm(particle_i.position_ - particle_j.position_),3));

    return particle_i.charge_ * E; 
}

vec PenningTrap::total_force_external(int i){

    Particle particle_i = particles.at(i);

    vec a = particle_i.charge_*external_E_field(particle_i.position_);
    vec b = particle_i.charge_*particle_i.velocity_;
    vec c = external_B_field(particle_i.position_);
    vec force = a + cross(b,c);
    return force; 
}

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

vec PenningTrap::total_force(int i){

    vec force = total_force_external(i) + total_force_particles(i);
    return force; 
}

void PenningTrap::evolve_forward_Euler(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){

        vec a = total_force(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}
void PenningTrap::evolve_Euler_Cromer(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){

        vec a = total_force(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}
void PenningTrap::evolve_Euler_Cromer_without_particle_interaction(double dt){

   int n = particles.size();

   for (int i=0; i < n; i++){
    //    Particle p = particles.at(i);

        vec a = total_force_external(i)/particles.at(i).mass_;
        particles.at(i).velocity_ = particles.at(i).velocity_ + a*dt;
        particles.at(i).position_ = particles.at(i).position_ + particles.at(i).velocity_*dt; 

   }

}
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