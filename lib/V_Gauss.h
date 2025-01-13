#ifndef V_GAUSS_H
#define V_GAUSS_H

#include "potential.h"

class V_Gauss : public LocalPotential {
  public:
    V_Gauss(const System& sys, double V0, double R) : LocalPotential(sys, "V_Gauss"), V0(V0), R(R), R2(R*R) {}
    
    double operator()(double r) const override {
      return V0 * std::exp(-r * r / R2);
    }

    double get(int ell, double p, double q) const override {
      if (ell == 0){
        double A = p * p + q * q;
        double B = 2 * p * q;
        return V0 * 4.0 * std::pow(M_PI, 1.5) * (R / B) 
          * std::exp(-0.25 * A * R2) * std::sinh(0.25 * B * R2);
      }
      else {
        return LocalPotential::get(ell, p, q);
      }
    }

  private:
    double V0;
    double R;
    double R2;
}; 

#endif // V_GAUSS_H
