#include "V_EFT_LO.h"
#include <map>

int main() {
    System sys(1.0, 1.0);
    V_EFT_LO potential(sys, 1.0, 2.0, 10000.);

    // Test p-space plotting
    potential.show("p", 0, {{"max", 10.0}, {"num", 50}});

    return 0;
}
