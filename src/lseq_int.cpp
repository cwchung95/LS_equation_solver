#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "system.h"
#include "potential.h"
#include "V_Gauss.h"
#include "V_Zero.h"
#include "V_EFT_LO.h"
#include "mesh.h"
#include "gaulegmesh.h"
#include "expgaulegmesh.h"
#include "operators.h"
#include "plot.h"

struct IntegrandParams {
  System* sys;
  LocalPotential* V;
  G0* g0;
  double k;
};

std::complex<double> Tmatrix_From_Integral(System& sys, double k, double mesh_min, double mesh_max, LocalPotential& V, G0& g0);
int get_shell_width();
void print_text_box(const std::vector<std::string>& lines);

std::string realToStringWithSigFigs(double x, int sigFigs); 
std::string complexToStringWithSigFigs(const std::complex<double>& c, int sigFigs);
std::vector<std::string> format_config(const std::map<std::string, std::string>& config, const std::vector<std::string>& config_order, const std::string& title, const std::string& filename);
std::string trim(const std::string& str);
std::map<std::string, std::string> parse_datacard(const std::string& filename);

double integrand(double p, void* params);
double integrand_real(double p, void* params); 
double integrand_imag(double p, void* params);  
double phase_shift_from_T(const System& sys, double k, std::complex<double> T);

int main(int argc, char* argv[]) {
  // Default configurations
  std::map<std::string, std::string> config = {
    {"mass", "1.0"},
    {"scale", "1.0"},
    {"mesh_scheme", "gauleg"},
    {"mesh_points", "16"},
    {"mesh_min", "0.0"},
    {"mesh_max", "4.0"},
    {"potential", "V_Gauss"},
    {"potential_params", "-4.0, 2.0"},
    {"k", "1.0"}
  };
  std::vector<std::string> config_order = {
    "mass", "scale", "mesh_scheme", "mesh_points", "mesh_min", "mesh_max", "potential", "potential_params", "k"
  };
  bool flag_debug = false;
  std::string datacard_file = "Default Setting - No datacard file provided";

  // Override the default configurations with the datacard file
  if (argc > 3 || (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
    std::cerr << "Usage: " << argv[0] << " [datacard_file] (-v/--verbose)" << std::endl;
    return 1;
  }
  else if (argc == 2) {
    datacard_file = argv[1];
    try{
      auto datacard_config = parse_datacard(datacard_file);
      for (const auto& [key, value] : datacard_config) {
        config[key] = value;
      }
    } 
    catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      return 1;
    }
  }
  else if (argc == 3 && (std::string(argv[2]) == "-v" || std::string(argv[2]) == "--verbose")) {
    flag_debug = true;
    datacard_file = argv[1];
    try {
      auto datacard_config = parse_datacard(datacard_file);
      for (const auto& [key, value] : datacard_config) {
        config[key] = value;
      }
    } 
    catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      return 1;
    }
  }

  // Extract system parameters
  double scale = std::stod(config["scale"]);
  double mass = std::stod(config["mass"]);

  //Extract mesh parameters
  std::string mesh_scheme = config["mesh_scheme"];
  size_t mesh_points = std::stoi(config["mesh_points"]);
  double mesh_min = std::stod(config["mesh_min"]);
  double mesh_max = std::stod(config["mesh_max"]);

  // Extract potential parameters
  std::string potential_type = config["potential"];
  std::vector<double> potential_params;
  std::istringstream params_stream(config["potential_params"]);
  std::string param;
  while (std::getline(params_stream, param, ',')) {
    potential_params.push_back(std::stod(trim(param)));
  }

  // Extract k
  double k = std::stod(config["k"]);

  // Initialize the system
  System sys;
  sys.setScale(scale);
  sys.setMass(mass);
  // std::cout << sys << std::endl;

  // Create the mesh
  Mesh* mesh = nullptr;
  if (mesh_scheme == "gauleg") {
    mesh = new GaulegMesh(mesh_points, mesh_min, mesh_max);
  } 
  else if (mesh_scheme == "expgauleg") {
    mesh = new ExpGaulegMesh(mesh_points, mesh_min, mesh_max);
  }
  else {
    std::cerr << "Error: Invalid mesh scheme: " << mesh_scheme << std::endl;
    return 1;
  }

  if (!mesh) {
    std::cerr << "Error: Mesh initialization failed!" << std::endl;
    return 1;
  }

  // Create the potential
  LocalPotential* V = nullptr;
  if (potential_type == "V_Gauss") {
    if (potential_params.size() != 2) {
      throw std::invalid_argument("Error: Invalid number of potential parameters for V_Gauss!");
    }
    V = new V_Gauss(sys, potential_params[0], potential_params[1]);
  }
  else if (potential_type == "V_Zero") {
    if (potential_params.size() != 0) {
      throw std::invalid_argument("Error: Invalid number of potential parameters for V_Zero!");
    }
    V = new V_Zero(sys);
  }
  else if (potential_type == "V_EFT_LO_unnat") {
    if (potential_params.size() != 3) {
      throw std::invalid_argument("Error: Invalid number of potential parameters for V_EFT_LO_unnat!");
    }
    V = new V_EFT_LO_unnat(sys, potential_params[0], potential_params[1], potential_params[2]);
  }
  else if (potential_type == "V_EFT_LO_nat") {
    if (potential_params.size() != 3) {
      throw std::invalid_argument("Error: Invalid number of potential parameters for V_EFT_LO_nat!");
    }
    V = new V_EFT_LO_nat(sys, potential_params[0], potential_params[1], potential_params[2]);
  }
  else {
    throw std::invalid_argument("Error: Invalid potential type: " + potential_type);
  }

  // Show the potential
  if (flag_debug) {
    V->show("r", 0, {{"max", mesh_max}});
    V->show("p", 0, {{"max", mesh_max}});
  }

  // Display configuration
  print_text_box(format_config(config, config_order, "CONFIGURATION", datacard_file)); 

  // Define the Green's function
  G0 G0(sys);

  // Calculate the T-matrix
  std::complex<double> T = Tmatrix_From_Integral(sys, k, mesh_min, mesh_max, *V, G0);
  std::cout << "T-matrix: " << complexToStringWithSigFigs(T, 4) << std::endl;

  // Calculate the phase shift
  double delta = phase_shift_from_T(sys, k, T);
  std::cout << "Phase shift: " << realToStringWithSigFigs(delta, 4) << " degrees" << std::endl;

  return 0;
}

std::complex<double> Tmatrix_From_Integral(System& sys, double k, double mesh_min, double mesh_max, LocalPotential& V, G0& g0) {
  // Define integrand
  IntegrandParams ip;
  ip.V   = &V;
  ip.g0  = &g0;
  ip.sys = &sys;
  ip.k   = k;

  // Set up the GSL integration workspace
  const size_t max_subdivisions = 10000;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(max_subdivisions);

  // Set up the GSL function
  gsl_function fr, fi;
  fr.function = integrand_real;
  fi.function = integrand_imag;
  fr.params = &ip;
  fi.params = &ip;

  // Perform the Integration
  double a(mesh_min), b(mesh_max), eps_abs(1e-6), eps_rel(1e-6);

  double resultRe, errorRe, resultIm, errorIm;

  int key=1;

  int statusRe = gsl_integration_qag(
    &fr, a, b, 
    eps_abs, eps_rel, 
    max_subdivisions, key,
    workspace, 
    &resultRe, &errorRe
  );

  int statusIm = gsl_integration_qag(
    &fi, a, b, 
    eps_abs, eps_rel, 
    max_subdivisions, key,
    workspace,  
    &resultIm, &errorIm
  );

  std::complex<double> result(resultRe, resultIm);

  // print the results with error
  std::cout << "Real part: " << resultRe << " +/- " << errorRe << std::endl;
  std::cout << "Imaginary part: " << resultIm << " +/- " << errorIm << std::endl;

  // Free the workspace
  gsl_integration_workspace_free(workspace);

  // Calculate the T-matrix
  std::complex<double> T = 4 * M_PI / sys.getMass() * 1.0 / (V.get(0, k, k) - result);
  return T;
}

double integrand(double p, void* params) {
  IntegrandParams* ip = static_cast<IntegrandParams*> (params);
  return (p * p) * ip->V->get(0, ip->k, p) * (*ip->g0)(ip->sys->e_from_k(ip->k) , p)/ (2 * M_PI);
};

double integrand_imag(double p, void* params) {
  IntegrandParams* ip = static_cast<IntegrandParams*> (params);
  double eps = 1e-8;
  double A = ip->sys->getMu()*ip->sys->e_from_k(ip->k)-p*p;
  double denom = A*A + eps*eps;

  // Im[G0] = -eps /denom;
  double G0_imag = - eps / denom;

  return p * p * ip->V->get(0, ip->k, p) * ip->sys->getMu() * G0_imag / (2 * M_PI);
};

double integrand_real(double p, void* params) { 
  IntegrandParams* ip = static_cast<IntegrandParams*> (params);
  double eps = 1e-8;
  double A = ip->sys->getMu()*ip->sys->e_from_k(ip->k)-p*p;
  double denom = A*A + eps*eps;

  // Re[G0] = A /denom;
  double G0_real = A / denom;

  return p * p * ip->V->get(0, ip->k, p) * ip->sys->getMu() * G0_real / (2 * M_PI);
};

// Function to convert a real number to a string with a specified number of significant figures
std::string realToStringWithSigFigs(double x, int sigFigs) {
    std::ostringstream oss;
    oss << std::setprecision(sigFigs) << x;
    return oss.str();
}

// Function to convert a complex number to a string with a specified number of significant figures
std::string complexToStringWithSigFigs(const std::complex<double>& c, int sigFigs) {
    std::ostringstream oss;
    oss << std::setprecision(sigFigs);
    
    if (c.real() >= 0) {
      oss << "(+";
    } else {
      oss << "(-";
    }

    oss << std::abs(c.real()) << ")+";
    if (c.imag() >= 0) {
        oss << "(+";
    } else {
        oss << "(-";
    }
    oss << std::abs(c.imag()) << ")i";
    
    return oss.str();
}

// Function to trim whitespace from a string
std::string trim(const std::string& str) {
  size_t first = str.find_first_not_of(" \t");
  size_t last = str.find_last_not_of(" \t");
  return (first == std::string::npos) ? "" : str.substr(first, (last - first + 1));
}

// Function to get the shell width
int get_shell_width() {
  FILE* pipe = popen("tput cols", "r");
  if (!pipe) {
    std::cerr << "Warning: Unable to get the shell width!" << std::endl;
    return 80;
  }

  char buffer[128];
  if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
    pclose(pipe);
    try {
      int width = std::stoi(buffer);
      if (width > 0) {
        return width;
      }
    }
    catch (const std::exception& e) {
      std::cerr << "Error: Invalid output from tput cols: " << buffer << std::endl;
    }
  } 
  else {
    pclose(pipe);
  }

  return 80;
}

// Function to print a dynamic text box
void print_text_box(const std::vector<std::string>& lines) {
  int width = get_shell_width();
  int box_width = width - 4;

  // Ensure the shell width is valid
  if (width <= 50) width = 60;

  // Top border
  std::cout << "+" << std::string(width - 2, '-') << "+" << std::endl;
  std::cout << "|" << std::string(width - 2, ' ') << "|" << std::endl;

  // Content lines
  for (const auto& line : lines) {
    int padding = (box_width - line.size()) / 2;
    std::cout << "| " << std::string(std::max<std::size_t>(padding, 0), ' ') << line
      << std::string(std::max<std::size_t>(box_width - padding - line.size(), 0), ' ') << " |" << std::endl;
  }

  // Bottom border
  std::cout << "|" << std::string(width - 2, ' ') << "|" << std::endl;
  std::cout << "+" << std::string(width - 2, '-') << "+" << std::endl;
}

// Function to format the configuration
std::vector<std::string> format_config(const std::map<std::string, std::string>& config, const std::vector<std::string>& config_order, const std::string& title, const std::string& filename) {
  std::vector<std::string> config_lines;

  int title_width = title.size() + 2;

  std::ostringstream oss;
  oss << "+" << std::string(title_width, '-') << "+";
  config_lines.push_back(oss.str());
  oss.str("");
  oss << "| "<< title << " |";
  config_lines.push_back(oss.str());
  oss.str("");
  oss << "+" << std::string(title_width, '-') << "+";
  config_lines.push_back(oss.str());
  oss.str("");
  
  int key_width = 20;

  for (const auto& key : config_order) {
    std::string value = config.at(key);
    int value_width = 55 - key_width - 3;
    std::ostringstream oss;

    oss << std::setw(key_width) << std::left << key << ": " << std::right << std::setw(value_width) << value;
    config_lines.push_back(oss.str());
  }

  if (title == "CONFIGURATION") {
    config_lines.push_back("");
    config_lines.push_back("Datacard: " + filename);
    config_lines.push_back("");
  }
  return config_lines;
}

// Function to parse a datacard file
std::map<std::string, std::string> parse_datacard(const std::string& filename) {
  std::map<std::string, std::string> config;
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Error: Unable to open the file " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    // Ignore comments and empty lines
    line = trim(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // Split the line into key and value
    size_t pos = line.find("=");
    if (pos == std::string::npos) {
      throw std::runtime_error("Error: Invalid line in the datacard file: " + line);
    }

    std::string key = trim(line.substr(0, pos));
    std::string value = trim(line.substr(pos + 1));
    config[key] = value;
  }
  
  return config;
}


std::vector< std::complex<double> > real_to_complex(const std::vector<double>& real_vec) {
  std::vector< std::complex<double> > complex_vec;
  for (double x : real_vec) {
    complex_vec.push_back(std::complex<double>(x, 0.0));
  }
  return complex_vec;
}

// Get the phase shift
double phase_shift_from_T(const System& sys, double k, std::complex<double> T) {
  std::complex<double> cot_delta = -2.0 * M_PI / (k * sys.getMu()) * (1.0 / T) + std::complex<double> (0,1);
  double delta = std::atan(1.0 / cot_delta.real());
  //delta = std::atan(T.imag()/T.real());
  return delta * 180.0 / M_PI;
}
