#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <cmath>

class System {
  public:
    System(double mass = 1.0, double scale = 1.0);

    double e_from_k(double k) const;
    double k_from_e(double e) const;
    double getMass() const {
      return mass;
    }

    std::string repr() const;

  private:
    double mass;
    double scale;
    double mu;
};


#endif // SYSTEM_H
