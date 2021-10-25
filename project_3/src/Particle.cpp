#include "Particle.hpp"

// Class for modeling a particle with its properties.
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
