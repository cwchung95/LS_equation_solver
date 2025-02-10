#include "system.h"


System::System(double mass, double scale) : mass(mass), scale(scale) {
  //mu = mass / 2.0;
  mu = mass / 2.0;
}

double System::e_from_k(double k) const {
  return 0.5 * k * k / mu;
}

double System::k_from_e(double e) const {
  return std::sqrt(2.0 * mu * e);
}

std::string System::repr() const {
  return "System(mass=" + std::to_string(mass) + ", scale=" + std::to_string(scale) + ")";
}
