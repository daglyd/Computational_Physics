#include "Particle.hpp"

/// <summary>Particle object for modelling a particle</summary>
/// <param name="charge">double: Charge of the particle</param>  
/// <param name="mass">double: mass of the particle</param>  
/// <param name="positions">arma::vec: initial position of the particle (x,y,z)</param> 
/// <param name="velocity">arma::vec: initial velocities of the particle (vx,vy,vz)</param> 
Particle::Particle(double charge, double mass, vec position, vec velocity)
{
  
  charge_ = charge;
  mass_ = mass;
  position_ = position;
  velocity_ = velocity;
}

// Definion of function
void Particle::some_function()
{
  cout << "This is some function" << endl;
}
