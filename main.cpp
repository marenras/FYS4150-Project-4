#include "ising.h"


int main(int argc, char* argv[])
{
  int mc_cycles = 100000;
  bool random_lattice = false;

  Ising model(2, "2D_lattice_100000.dat", 1, 4, 0.2, mc_cycles, random_lattice);

  //Ising model(20, "test.dat");
  //model.InitializeLattice(1, random_lattice);
  //model.MonteCarloSample(mc_cycles);
}
