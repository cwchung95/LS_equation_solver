# LS equation solver

## How to build
You can compile the entire code with `compile.sh` script.
```bash
./compile.sh clean
./compile.sh
```
where `./compile.sh` will compile the code, and `./compile.sh clean` will clean up the executable and object files.

Also, if you cannot build the code with `./compile.sh`, you can build the code manually with the following commands:
```bash
mkdir build
cd build
cmake ..
make 
```
This will allow you to get the executables in `run` directory.

## Source code explanation
The source code can be found at `src` directory. Currently, there is only one source code: `lseq.cpp` which solves the Lippmann-Schwinger equation with given potential.
After build the code, you can find the executable in `run` directory, and run it with following command:
```bash
./run/lseq [datacard] (-v/--verbose)
```
where `datacard` is the input file that contains the parameters for the calculation, and `-v` or `--verbose` is the option for verbose output.
By modifying the datacard, you can change the parameters for the calculation. Follwing is the list of parameters you can modify in the datacard:

| Parameter name     | Type          | Explanation                                                                                                      | Default Value |
|--------------------|---------------|------------------------------------------------------------------------------------------------------------------|---------------|
| `mass`             | `double`      | mass of nucleon (in MeV)                                                                                         | `1.0`         |
| `scale`            | `double`      | scale of the model you want to calculate (hbar c)                                                                | `1.0`         |
| `mesh_scheme`      | `std::string` | scheme of mesh structures  (possible options: `gauleg`, `expgauleg`)                                             | `gauleg`      |
| `mesh_points`      | `int`         | number of discrete points between momentum range `mesh_min` and `mesh_max`                                       | `16`          |
| `mesh_min`         | `double`      | minimum momentum (in MeV)                                                                                        | `0.0`         |
| `mesh_max`         | `double`      | maximum momentum (in MeV)                                                                                        | `4.0`         |
| `potential`        | `std::string` | set forms of potential (in MeV)  (possible options: `V_Gauss`, `V_Zero`, `V_EFT_LO`)                             | `V_Gauss`     |
| `potential_params` | `std::string` | set potential's parameters, number of parameters depends on potential (`V_Gauss`: 2, `V_Zero`: 0, `V_EFT_LO`: 3) | `-4.0, 2.0`   |
| `k`                | `double`      | value of `k`                                                                                                     | `1.0`         |

Currently, the `V_EFT_LO` potential is not fully implemented, so it is recommended to use `V_Gauss` or `V_Zero` potential.

