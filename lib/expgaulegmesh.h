#ifndef EXPGAULEGMESH_H
#define EXPGAULEGMESH_H

#include <iostream>
#include "gaulegmesh.h"
#include <cmath>

class ExpGaulegMesh : public GaulegMesh {
  public:
    ExpGaulegMesh(int n, double p_min, double p_max) : GaulegMesh(n, std::log(1.0+p_min), std::log(1.0+p_max)) {
      for (int i = 0; i < n; i++) {
        points[i] = std::exp(points[i]) - 1.0;
        weights[i] *= (points[i] + 1.0);
      }
    }
};

#endif // EXPGAULEGMESH_H
