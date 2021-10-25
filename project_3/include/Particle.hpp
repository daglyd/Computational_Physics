#ifndef __particle__hpp__ 
#define __particle__hpp__ 

#include <armadillo>

using namespace std;
using namespace arma;

  

class Particle 
{
private:
  // Declaration of variables only accessible from within the class

public:
    
    double charge_, mass_;
    vec position_, velocity_;


    // Declaration of constructors, e.g.
    Particle(double charge, double mass, vec position, vec velocity);

    // Declaration of destructors, copy constructors, ...

    // Declarations of other class methods, e.g.
    void some_function();

}; // <-- Note that class bodies end with a semicolon!

#endif