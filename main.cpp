#include "ising.h"


int main(int argc, char* argv[])
{
  int mc_cycles = 100000;
  Ising model(2, "2D_lattice_100000.dat", 1, 4, 0.2, mc_cycles);
}
