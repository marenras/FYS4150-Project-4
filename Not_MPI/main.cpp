#include "ising.h"


int main(int argc, char* argv[])
{
    if (argc != 5 && argc != 7){
      std::cout << "Bad Usage: " << argv[0] <<
        " Please provide lattice dimension, MC cycles, initial temperature, final temperature, tempurate step and output filename. " <<
        "If you want to run for only one temperature, provide lattice dimension, MC cycles, temperature and output filename." << std::endl;
      exit(1);
    }

    std::string filename;
    int dim_lattice, MC_cycles;
    double T_min, T_max, dT, T;
    bool random_lattice = true;


    if (argc == 7) {
      dim_lattice = atoi(argv[1]);
      MC_cycles   = atoi(argv[2]);
      T_min       = atof(argv[3]);
      T_max       = atof(argv[4]);
      dT          = atof(argv[5]);
      filename    = argv[6];

      Ising model(dim_lattice, filename, T_min, T_max, dT, MC_cycles, random_lattice);
    }


    if (argc == 5) {
      dim_lattice = atoi(argv[1]);
      MC_cycles   = atoi(argv[2]);
      T           = atoi(argv[3]);
      filename    = argv[4];

      Ising model(dim_lattice, filename);
      model.InitializeLattice(T, random_lattice);
      model.MonteCarloSample(MC_cycles);
    }

}
