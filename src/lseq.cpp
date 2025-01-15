#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "system.h"
#include "potential.h"
#include "V_Gauss.h"
#include "mesh.h"
#include "gaulegmesh.h"
#include "operators.h"
#include "plot.h"

std::vector< std::vector<std::complex<double> > > kernel(const Mesh& mesh, const G0& G0, const Potential& V, double E);
std::vector< std::vector< std::complex<double> > > construct_I_minus_K(const size_t n, const std::vector< std::vector< std::complex<double> > > K);
std::vector< std::complex<double> > real_to_complex(const std::vector<double>& real_vec);
std::vector< std::complex<double> > linear_solve(const std::vector< std::vector< std::complex<double> > >& mat, 
                                                const std::vector< std::complex<double> >& vec);
std::vector< std::complex<double> > solve(const System& sys, Mesh& mesh, const G0& G0, const Potential& V, double k);
std::complex<double> solve_on_shell(const System& sys, Mesh& mesh, const G0& G0, const Potential& V, double k);
double phase_shift_from_T(const System& sys, double k, std::complex<double> T);

int main() {
  // Initialize the system
  System sys;
  std::cout << sys << std::endl;

  // Define the potential
  V_Gauss V(sys, -4.0, 2.0);

  // Show the potential
  V.show("r", 0, {{"max", 4.0}});
  V.show("p", 0, {{"max", 4.0}});

  // Create the mesh
  GaulegMesh mesh(128, 0.0, 4.0);

  // Define the Green's function
  G0 G0(sys);

  // Plot the Green's function
  double k = 1.0;
  std::vector<double> ps = mesh.ps();
  std::vector<double> G0_values;
  for (double p : ps) {
    G0_values.push_back(G0(sys.e_from_k(k), p));
  }

  // Solve the integral equation
  auto solution = solve(sys, mesh, G0, V, k);
  std::cout << "Solution: " << std::endl;
  for (std::complex<double> s : solution) {
    std::cout << s << " ";
  }

  // Calculate the T-matrix
  std::complex<double> T = solve_on_shell(sys, mesh, G0, V, k);
  std::cout << "\nT-matrix : " << T << std::endl;

  //Calculate the Phase Shift
  double delta = phase_shift_from_T(sys, k, T);
  std::cout << "Phase Shift: " << delta << " degrees" << std::endl;
}


std::vector< std::complex<double> > real_to_complex(const std::vector<double>& real_vec) {
  std::vector< std::complex<double> > complex_vec;
  for (double x : real_vec) {
    complex_vec.push_back(std::complex<double>(x, 0.0));
  }
  return complex_vec;
}

std::vector< std::complex<double> > linear_solve(const std::vector< std::vector< std::complex<double> > >& mat, 
                                                const std::vector< std::complex<double> >& vec) {
  size_t n = mat.size();

  // Check if the matrix is square and match with the vector size
  if (n == 0 || mat[0].size() != n || vec.size() != n) {
    throw std::invalid_argument("Matrix and vector sizes do not match!");
  }

  // Allocate memory for the GSL matrix and vector
  gsl_matrix_complex* gsl_mat = gsl_matrix_complex_alloc(n, n);
  gsl_vector_complex* gsl_vec = gsl_vector_complex_alloc(n);
  gsl_vector_complex* gsl_sol = gsl_vector_complex_alloc(n);
  gsl_permutation* perm = gsl_permutation_alloc(n);
  int signum;

  // Copy the matrix and vector to the GSL structures
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      gsl_complex z = gsl_complex_rect(mat[i][j].real(), mat[i][j].imag());
      gsl_matrix_complex_set(gsl_mat, i, j, z);
    }
  }

  for (size_t i = 0; i < n; i++) {
    gsl_complex z = gsl_complex_rect(vec[i].real(), vec[i].imag());
    gsl_vector_complex_set(gsl_vec, i, z);
  }

  // Solve the linear system
  gsl_linalg_complex_LU_decomp(gsl_mat, perm, &signum);
  gsl_linalg_complex_LU_solve(gsl_mat, perm, gsl_vec, gsl_sol);

  // Extract the solution
  std::vector< std::complex<double> > sol(n);
  for (size_t i = 0; i < n; i++) {
    gsl_complex z = gsl_vector_complex_get(gsl_sol, i);
    sol[i] = std::complex<double>(GSL_REAL(z), GSL_IMAG(z));
  }

  // Free the memory
  gsl_matrix_complex_free(gsl_mat);
  gsl_vector_complex_free(gsl_vec);
  gsl_vector_complex_free(gsl_sol);
  gsl_permutation_free(perm);

  return sol;
}

// Kernel Matrix
std::vector< std::vector<std::complex<double> > > kernel(const Mesh& mesh, const G0& G0, const Potential& V, double E) {
  double factor = 0.5 / (M_PI * M_PI);
  std::vector< std::vector<std::complex<double> > > K;

  for (size_t i=0; i < mesh.size(); i++ ) {
    std::vector<std::complex<double> > row;
    for (size_t j=0; j < mesh.size(); j++) {
      std::complex<double> value = factor * mesh.w(j) * mesh.p(j) * mesh.p(j) 
                                  * (mesh.is_pv(j) ?  G0.residue(mesh.p(j)) : G0(E, mesh.p(j)))
                                  * V.get(0, mesh.p(i), mesh.p(j));
      row.push_back(value);
    }
    K.push_back(row);
  }

  return K;
}

// Generate the T-matrix from K-matrix
std::vector< std::vector< std::complex<double> > > construct_I_minus_K(const size_t n, const std::vector< std::vector< std::complex<double> > > K) {
  std::vector< std::vector< std::complex<double> > > I_minus_K(n, std::vector< std::complex<double> >(n, 0.0));
  for (size_t i = 0; i < n; ++i) {
    for(size_t j = 0; j < n; ++j) {
      if(i==j) {
        I_minus_K[i][j] = 1.0 - K[i][j];
      } else {
        I_minus_K[i][j] = -K[i][j];
      }
    }
  }
  return I_minus_K;
}

// Solve the integral equation
std::vector< std::complex<double> > solve(const System& sys, Mesh& mesh, const G0& G0, const Potential& V, double k) {
  mesh.push_pv(k); 

  double E = sys.e_from_k(k);
  auto K = kernel(mesh, G0, V, E);
  std::vector<double> vec;

  size_t n = mesh.size();
  auto I_minus_K = construct_I_minus_K(n, K);

  for(double p : mesh.ps()) {
    vec.push_back(V.get(0, k, p));
  }
  std::vector< std::complex<double> > solution = linear_solve(I_minus_K, real_to_complex(vec));

  mesh.pop_pv();
  return solution;
}

// Solve on-shell
std::complex<double> solve_on_shell(const System& sys, Mesh& mesh, const G0& G0, const Potential& V, double k){
  int i0 = mesh.push_pv(k);

  double E = sys.e_from_k(k);
  auto K = kernel(mesh, G0, V, E);
  std::vector<double> vec;

  size_t n = mesh.size();
  auto I_minus_K = construct_I_minus_K(n, K);

  for (double p : mesh.ps()) {
    vec.push_back(V.get(0, k, p));
  }

  std::vector< std::complex<double> > solution = linear_solve(I_minus_K, real_to_complex(vec));

  mesh.pop_pv();
  return solution[i0];
}

// Get the phase shift
double phase_shift_from_T(const System& sys, double k, std::complex<double> T) {
  std::complex<double> cot_delta = -2.0 * M_PI / (k * sys.getMu()) * (1.0 / T) + std::complex<double> (0,1);
  double delta = std::atan(1.0 / cot_delta.real());
  return delta * 180.0 / M_PI;
}
