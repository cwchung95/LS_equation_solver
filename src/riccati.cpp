#include <gsl/gsl_sf_bessel.h>
#include <cmath>


#include "riccati.h"

double riccati_j(int ell, double z) {
  return z * gsl_sf_bessel_jl(ell, z);
}

double riccati_n(int ell, double z) {
  return -z * gsl_sf_bessel_yl(ell, z);
}
