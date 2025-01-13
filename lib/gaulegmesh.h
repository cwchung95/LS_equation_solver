#ifndef GAULEGMESH_H
#define GAULEGMESH_H

#include "mesh.h"
#include <gsl/gsl_integration.h>
#include <cmath>

class GaulegMesh : public Mesh {
  public: 
    GaulegMesh(int n, double p_min, double p_max) : Mesh(n, p_min, p_max) {
      points.resize(n);
      weights.resize(n);
      legendre_gauss(n, p_min, p_max, points, weights);
    }

  private:
    void legendre_gauss(int n, double p_min, double p_max, std::vector<double>& points, std::vector<double>& weights) { 
      gsl_integration_glfixed_table* table = gsl_integration_glfixed_table_alloc(n);

      for (int i = 0; i < n; ++i) {
        double xi, wi;

        gsl_integration_glfixed_point(p_min, p_max, i, &xi, &wi, table);

        points[i] = xi;
        weights[i] = wi;
      }

      gsl_integration_glfixed_table_free(table);
    }
};

#endif // GAULEGMESH_H
