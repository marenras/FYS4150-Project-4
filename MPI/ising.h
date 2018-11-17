#include <armadillo>
#include <string>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>


class Ising
{
public:
  // Constructor
  Ising(int dim_lattice);

  // Initialize variables
  int dim_lattice;
  int N_spins;
  int MC_cycles;
  int accepted_states;
  double temperature;
  double T_min;
  double T_max;
  double dT;
  bool temp_header;
  bool random_lattice;

  arma::mat lattice;
  arma::vec expectation_values;
  arma::vec energy_difference;
  std::ofstream ofile;

  double energy;
  double mean_energy;
  double energy_variance;
  double specific_heat;

  double magnetization;
  double mean_magnetization;
  double susceptibility;
  double mean_absolute_magnetization;

  // Initialize random number generator
  //std::random_device rd; seed?
  std::mt19937 generator;
  std::uniform_real_distribution<double> RNG;


  // Initialize functions
  void InitializeLattice(double temperature, bool random_lattice);
  void Metropolis();
  int PBC(int i, int limit, int add);
};
