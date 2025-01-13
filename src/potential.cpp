#include "potential.h"
#include "plot.h"
#include "riccati.h"
#include <cmath>
#include <stdexcept>
#include <functional>

double integrate(std::function<double(double)> func, double a, double b) {
  int n = 1e6;
  double step = (b - a) / n;
  double result = 0.0;
  double x1 = a;
  double x2 = a + step;
  for (int i = 0; i < n; i++) {
    x1 = x1 + step;
    x2 = x2 + step;
    result += 0.5 * (func(x1) + func(x2)) * step; 
  }
  return result;
}

double LocalPotential::get(int ell, double p, double q) const {
  auto func = [&](double r){
    return riccati_j(ell, p*r) * (*this)(r) * riccati_j(ell, q*r);
  };
  return 4.0 * M_PI / (q * p) * integrate(func, 0.0, 1e6);
}

void LocalPotential::show(const std::string& rep, int ell, const std::map<std::string, double>& options) const {
  Plot plot;

  if(rep == "r"){
    double max = options.count("max") ? options.at("max") : 10.0;
    int num = options.count("num") ? static_cast<int>(options.at("num")) : 32;

    std::vector<double> r_values, v_values;
    double r = (max / num);
    for (int i = 0; i < num; i++) {
      r_values.push_back(r);
      v_values.push_back((*this)(r));
      r = r + max / num;
    }

    plot.plot1D(r_values, v_values, "Potential in r-space");
  }
  else if (rep == "p") {
    double max = options.count("max") ? options.at("max") : 10.0;
    int num = options.count("num") ? static_cast<int>(options.at("num")) : 32;

    std::vector<double> px_values, py_values, v_values;

    double p = max / num;
    double q = max / num;
    for (int i = 0; i < num; i++) {
      for (int j = 0; j < num; j++) {
        px_values.push_back(p);
        py_values.push_back(q);
        v_values.push_back(this->get(ell, p, q));

        q = q + max / num;
      }
      q = max / num;
      p = p + max / num;
    }
    plot.plot2D(px_values, py_values, v_values, "Potential in p-space");
    //throw std::runtime_error("p-space plotting not yet implemented!");
  }
  else {
    throw std::runtime_error("Unsupported potential representation!");
  }
}
