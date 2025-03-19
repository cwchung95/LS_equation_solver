#ifndef OPERATORS_H
#define OPERATORS_H

#include <iostream>
#include <string>
#include <vector>

#include "system.h"

class Operators {
  public:
    Operators(const System& sys) {}

    virtual ~Operators() {}

  private:
    double mass;
    double mu;
};

class G0 : public Operators {
  public:
    G0(const System& sys) : Operators(sys), sys(sys), mass(sys.getMass()), mu(sys.getMu()) {}

    double operator()(double E, double q) const {
      return mu / (mu * E - q * q);
    }
    
    double residue(double q0) const {
      return -mu / (q0);
    }

  private:
    double mass;
    double mu;
    const System& sys;
};

#endif // OPERATORS_H
