#ifndef V_ZERO_H
#define V_ZERO_H

#include "potential.h"

class V_Zero : public LocalPotential {
  public:
    V_Zero(const System& sys) : LocalPotential(sys, "V_Zero") {}

    double operator()(double r) const override {
      return 0.0;
    }

};

#endif // V_ZERO_H
