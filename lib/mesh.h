#ifndef MESH_H
#define MESH_H

#include <vector>
#include <complex>
#include <cmath>
#include <numeric>

class Mesh {
  protected:
    int n;
    double p_min, p_max;
    std::vector<double> points, weights, pv_points;
    std::vector< std::complex<double> > pv_weights;

  public:
    Mesh(const int n, const double p_min, const double p_max) : n(n), p_min(p_min), p_max(p_max) {}

    virtual ~Mesh() {}

    std::vector<double> ps() const{
      std::vector<double> result = points;
      result.insert(result.end(), pv_points.begin(), pv_points.end());
      return result;
    }

    std::vector<std::complex<double> > ws() const {
      std::vector<std::complex<double> > result(weights.begin(), weights.end());
      result.insert(result.end(), pv_weights.begin(), pv_weights.end());
      return result;
    }

    double p(int i) const { return ps()[i]; }
    std::complex<double> w(int i) const { return ws()[i]; }

    size_t size() const { return n + pv_points.size(); }

    int n_pv(int i=0) const { return n + i; }

    int push_pv(double p0, bool ipi = true) {
      std::complex<double> R = 0.0;

      for (size_t i = 0; i < size(); ++i) {
        double p_i = p(i);
        R += w(i) / (p0 - p_i);
      }
      R += std::log((p_max - p0) / (p0 - p_min));
      if (ipi) {
        R += std::complex<double>(0.0, M_PI);
      }

      pv_points.push_back(p0);
      pv_weights.push_back(R);
      
      return size() - 1;

    }

    void pop_pv() {
      if (!pv_points.empty()) {
        pv_points.pop_back();
        pv_weights.pop_back();
      }
    }
};

#endif // MESH_H
