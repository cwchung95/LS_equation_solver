#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <cmath>

class System {
  public:
    System(double mass = 1.0, double scale = 1.0);

    double e_from_k(double k) const;
    double k_from_e(double e) const;

    double setMass(double m) {
      mass = m;
      return mass;
    }
    double getMass() const {
      return mass;
    }
    double setMu(double m) {
      mu = m;
      return mu;
    }
    double getMu() const {
      return mu;
    }
    double setScale(double s) {
      scale = s;
      return scale;
    }
    double getScale() const {
      return scale;
    }

    std::string repr() const;

    friend std::ostream& operator<<(std::ostream& os, const System& sys) {
      os << sys.repr();
      return os;
    }

  private:
    double mass;
    double scale;
    double mu;
};


#endif // SYSTEM_H
