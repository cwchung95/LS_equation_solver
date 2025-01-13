#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

class Plot {
  public:
    Plot () {
      gnuplotPipe = popen("gnuplot -persistent", "w");
      if(!gnuplotPipe){
        throw std::runtime_error("Could not open gnuplot pipe.");
      }
      command("set term qt font 'Arial,10'");
    }

    ~Plot() {
      if (gnuplotPipe) {
        pclose(gnuplotPipe);
      }
    }

    void command(const std::string& cmd) {
      fprintf(gnuplotPipe, "%s\n", cmd.c_str());
      fflush(gnuplotPipe);
    }
    
    void remove_temp_file() {
      std::string cmd = "rm temp_data.dat";
      system(cmd.c_str());
    }

      
    void plot1D(const std::vector<double>& x, const std::vector<double>& y, const std::string& title = "Plot") {
      if(x.size() != y.size()) {
        throw std::runtime_error("X and Y data sizes are not matching");
      }

      std::string tempFileName = "temp_data.dat";
      std::ofstream dataFile(tempFileName);
      for ( size_t i = 0; i < x.size(); ++i ) {
        dataFile << x[i] << " " << y[i] << "\n";
      }
      dataFile.close();

      std::string cmd = "plot '"+ tempFileName + "' with lines title '" + title + "'";
      command(cmd);

      //remove_temp_file();
    }

    void plot2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::string& title = "Plot") {
      if(x.size() != y.size()) {
        throw std::runtime_error("data sizes are not matching");
      }
      else if(x.size() != z.size()) {
        throw std::runtime_error("data sizes are not matching");
      }

      std::string tempFileName = "bin/temp_data.dat";
      std::ofstream dataFile(tempFileName);
      double y_max = *std::max_element(y.begin(), y.end());
      for ( size_t i = 0; i < x.size(); ++i ) { 
        if (std::isinf(z[i])||std::isnan(z[i])) continue;
        dataFile << x[i] << " " << y[i] << " " << z[i] << "\n";
        if (y[i] == y_max) {
          dataFile << "\n";
        }
      } 
      dataFile.close();

      std::string cmd = "splot '"+ tempFileName + "' using 1:2:3 with pm3d";
      
      command("set xlabel 'p'");
      command("set ylabel 'p\\''");
      command("set title 'Potential Representation'");
      command("set view map");
      command("unset key");
      command(cmd);

      //remove_temp_file();
    }


  private:
    FILE* gnuplotPipe;
};

#endif // PLOT_H
