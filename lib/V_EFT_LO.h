#ifndef V_EFT_LO_H
#define V_EFT_LO_H

#include "potential.h"
#include "system.h"

class V_EFT_LO : public LocalPotential {
  public:
    V_EFT_LO(const System& sys, double a, double mu, double lamd) : LocalPotential(sys, "V_EFT_LO"), a(a), mu(mu), lamd(lamd) {} 

    double operator()( double r ) const override {
      double a_MeV = sys.fm_to_MeV_inv(a);
      // std::cout << a_MeV << std::endl;
      return 32 * std::pow(M_PI, 2) * a_MeV / ( 1 - 2 * a_MeV * mu / M_PI ) * gaussian(r);
    }

    double gaussian(double r) const {
      return std::exp(-r * r / (2 * std::pow(mu, 2)))/(std::pow(2*M_PI, 1.5) * std::pow(mu, 3));
    }

    double get(int ell, double p, double q) const override {
      return LocalPotential::get(ell, p ,q);
      if (ell == 0){
        double A = p * p + q * q;
        double B = 2 * p * q;
        return V_EFT_LO::operator()(q) 
          * std::exp(-0.25 * A * std::pow(mu, 2)) * std::sinh(0.25 * B * std::pow(mu, 2));
      }
      else {
        return LocalPotential::get(ell, p, q);
      }
    }
    
  private:
    double a;
    double mu;
    double lamd;
}; 

#endif // V_EFT_LO_H
