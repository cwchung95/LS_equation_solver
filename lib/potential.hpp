#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <map>

#include "riccati.h"
#include "system.h"


class Potential {
  protected:
    System sys;
    std::string name;

  public:
    Potential(const System& sys, const std::string& name) : sys(sys), name(name) {}

    virtual ~Potential() {}
    virtual double get(int ell, double p, double q) const = 0;
    virtual void show(const std::string& rep = "p", int ell = 0) const {
      throw std::runtime_error("Abstract method call!");
    }

    friend std::ostream& operator<<(std::ostream& os, const Potential& potential) {
      os << potential.name;
      return os;
    }

};

class LocalPotential : public Potential {
  public:
    LocalPotential(const System& sys, const std::string& name)
      : Potential(sys, name) {}

    virtual double operator()(double r) const = 0;
    virtual double get(int ell, double p, double q) const override;
    virtual void show(const std::string& rep, int ell, const std::map<std::string, double>& options) const;
};

#endif // POTENTIAL_H
