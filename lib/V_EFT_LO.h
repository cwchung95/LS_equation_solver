#ifndef V_EFT_LO_H
#define V_EFT_LO_H

#include <gsl/gsl_sf_bessel.h>
#include "potential.h"
#include "system.h"

class V_EFT_LO_unnat : public LocalPotential {
  public: 
    // Constructor: get the system and the parameters (scattering length, effective range, an arbitrary energy scale) 
    V_EFT_LO_unnat( const System& sys, double a, double r0, double mu ) : LocalPotential( sys, "V_EFT_LO" ), a_MeV( sys.fm_to_MeV_inv(a) ), r0_MeV( sys.fm_to_MeV_inv(r0) ) {} 
    
    // Calculate the potential at a given distance r
    double operator()( double r ) const override {
      double C0 = 4 * M_PI / sys.getMass() * ( 1 / ( -mu + 1/a_MeV ) );
      return C0 * gaussian(r, r0_MeV);
    }

    // Gaussian function
    double gaussian( double r, double r0 ) const {
      return exp(-r*r/(r0*r0)) / (std::pow(2*M_PI, 1.5) * r0 * r0 * r0);
    }
    
    // Calculate V(p,q) from V(r) using the Lippman-Schwinger equation
    double get( int ell, double p, double q ) const override {
      double C0 = 4 * M_PI / sys.getMass() * ( 1 / ( - mu + 1/a_MeV ) );
      return C0;
      if ( ell == 0 ) {
        double A = p * p + q * q;
        double B = 2 * p * q;
        double ret = C0 * std::exp( - A * r0_MeV * r0_MeV ) * gsl_sf_bessel_i0_scaled( 2 * A * r0_MeV );
        return ret;
      }
      else {
        return LocalPotential::get(ell, p, q);
      }
    }

  private:
    double a_MeV, r0_MeV, mu;
};

class V_EFT_LO_nat : public LocalPotential {
  public:
    // Constructor: get the system and the parameters (scattering length, effective range, an arbitrary energy scale)
    V_EFT_LO_nat( const System& sys, double a, double r0, double mu ) : LocalPotential( sys, "V_EFT_LO" ), a_MeV( sys.fm_to_MeV_inv(a) ), r0_MeV( sys.fm_to_MeV_inv(r0) ) {}

    // Calculate the potential at a given distance r
    double operator()( double r ) const override {
      double C0 = 4 * M_PI * a_MeV / sys.getMass();
      return C0 * gaussian(r, r0_MeV);
    }

    // Gaussian function
    double gaussian( double r, double r0 ) const {
      return exp(-r*r/(r0*r0)) / (std::pow(2*M_PI, 1.5) * r0 * r0 * r0);
    }

    // Calculate V(p,q) from V(r) using the Lippman-Schwinger equation
    double get( int ell, double p, double q ) const override {
      double C0 = 4 * M_PI * a_MeV / sys.getMass();
      return C0;
      if ( ell == 0 ) {
        double A = p * p + q * q;
        double B = 2 * p * q;
        double ret = C0 * std::exp( - A * r0_MeV * r0_MeV ) * gsl_sf_bessel_i0_scaled( 2 * A * r0_MeV );
        return ret;
      }
      else {
        return LocalPotential::get(ell, p, q);
      }
    }
  
  private:
    double a_MeV, r0_MeV, mu;
};

#endif // V_EFT_LO_H

