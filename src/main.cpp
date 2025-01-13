#include "system.h"
#include "operators.h"
#include <iostream>

int main() {
  System sys(1.0, 1.0);

  G0 g0(sys);

  double E = 2.0;
  double q = 1.5;

  std::cout<< "G0 = " << g0(E, q) << std::endl;

  double q0 = 1.5;
  std::cout<< "Residue(q0) = " << g0.residue(q0) << std::endl;
}
