#include "ising.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <random>
 
Ising::Ising(int dim_lattice, std::string filename, bool temp_header)
{
  // Setting "self" values
  this->dim_lattice = dim_lattice;
  this->filename = filename;

  // Creating emty lattice
	N_spins = dim_lattice*dim_lattice;
	lattice = arma::mat(dim_lattice, dim_lattice);
	expectation_values = arma::vec(5);
	energy_difference = arma::vec(17);

	// Energy
	energy = 0;
	mean_energy = 0;
	energy_variance = 0;
	specific_heat = 0;

	// Magnetization
	magnetization = 0;
	mean_magnetization = 0;
	susceptibility = 0;
	mean_absolute_magnetization = 0;

  temp_header = false;
  WriteHeader(temp_header);
}

Ising::Ising(int dim_lattice, std::string filename, double T_min, double T_max, double dT, int mc_cycles)
{
  // Setting "self" values
  this->dim_lattice = dim_lattice;
  this->filename = filename;
  this->T_min = T_min;
  this->T_max = T_max;
  this->dT = dT;

  // Creating emty lattice
	N_spins = dim_lattice*dim_lattice;
	lattice = arma::mat(dim_lattice, dim_lattice);
	expectation_values = arma::vec(5);
	energy_difference = arma::vec(17);

	// Energy
	energy = 0;
	mean_energy = 0;
	energy_variance = 0;
	specific_heat = 0;

	// Magnetization
	magnetization = 0;
	mean_magnetization = 0;
	susceptibility = 0;
	mean_absolute_magnetization = 0;

  temp_header = true;
  WriteHeader(temp_header);

  for (double T = T_min; T <= T_max; T+=dT){
    InitializeLattice(T);
    MonteCarloSample(mc_cycles);
  }
}

void Ising::InitializeLattice(double temperature)
{
  // Setting "self" temperature
  this->temperature = temperature;


  // Reseting values
  energy = 0;
  magnetization = 0;
  accepted_states = 0;
  expectation_values.zeros();


  // Initialzing seed and the Mersienne random generator
	std::random_device rd;
	std::mt19937_64 generator(rd());
	std::uniform_real_distribution<double> RNG(0.0, 1.0);


  // Setup spin matrix with all spin values equal to 1 (up)
  // and initialize magnetization (sum of all spins)
  for (int x = 0; x < dim_lattice; x++) {
    for (int y = 0; y < dim_lattice; y++) {
      lattice(x,y) = 1.0;
      magnetization += lattice(x,y);
    }
  }

  // Setup initial energy
  for (int x = 0; x < dim_lattice; x++) {
    for (int y = 0; y < dim_lattice; y++) {
      // sum over nearest neigbours with periodic boundary conditions
      energy -= lattice(x,y) * (lattice(PBC(x, dim_lattice,-1), y)
      + lattice(x, PBC(y, dim_lattice,-1)));
    }
  }

  // Precomputing possible initial energies for a given temperature
	for (int i = -8; i <= 8; i++){
		energy_difference(i+8) = 0;
	}

	for (int i = -8; i <= 8; i += 4){
		energy_difference(i+8) = exp(-i/temperature);
	}
}



int Ising::PBC(int i, int limit, int add) {
  // Finds the nearest neighbours with periodic boundary conditions
  return (i + limit + add) % (limit);
}



void Ising::Metropolis()
{
	// The Metropolis algorithm.
	for (int i = 0; i < N_spins; i++)
		{
			// Removing potensial bias
			int rand_x = RNG(generator)*dim_lattice;
			int rand_y = RNG(generator)*dim_lattice;

			double delta_e =  2*lattice(rand_x, rand_y)*
        (lattice(rand_x, PBC(rand_y, dim_lattice, -1)) +
         lattice(PBC(rand_x, dim_lattice, -1), rand_y) +
         lattice(rand_x, PBC(rand_y, dim_lattice, 1)) +
         lattice(PBC(rand_x, dim_lattice, 1), rand_y));

			double rand_condition = RNG(generator);

			// Update variables if Metropolis condition is met
			if (rand_condition <= energy_difference(delta_e + 8))
			{
        // Flip one spin and accept new spin config
				lattice(rand_x, rand_y) *= -1;
				magnetization += 2 * lattice(rand_x, rand_y);
				energy += delta_e;
				accepted_states++;
			}
		}
}


void Ising::MonteCarloSample(int N)
{
	// Starts the Monte-Carlo sampling of the Ising-model for N MC-cycles using
	// the Metropolis algorithm.

	MC_cycles = N;

	// Starts the Monte-Carlo sampling
	for (int i = 1; i < MC_cycles; i++)
	{
		Metropolis();

		expectation_values(0) += energy;
		expectation_values(1) += energy * energy;
		expectation_values(2) += magnetization;
		expectation_values(3) += magnetization * magnetization;
		expectation_values(4) += fabs(magnetization);

    WriteToFile(i);
	}

}



void Ising::WriteHeader(bool temp_header)
{
  using namespace std;

  ofile.open("Data/"+filename, ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  if (temp_header == true)
   {
       ofile << setprecision(4) << T_min;
       ofile << setw(8) << setprecision(4) << T_max;
       ofile << setw(8) << setprecision(4) << dT;
   }
   ofile << "\n" << "\n";
}


void Ising::WriteToFile(int current_cycle)
{
  using namespace std;
  ofile.open("Data/"+filename, ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  // Normalizing expectation values
  double norm = 1.0 / ((double) (current_cycle));
  arma::vec normalized_expectation_values = expectation_values*norm;

  // Variance calculations
	energy_variance = (normalized_expectation_values(1) - normalized_expectation_values(0)*normalized_expectation_values(0)) / N_spins;
	susceptibility = (normalized_expectation_values(3) - normalized_expectation_values(2)*normalized_expectation_values(2)) / (N_spins * temperature);

  mean_energy = normalized_expectation_values(0) / N_spins;
  mean_magnetization = normalized_expectation_values(2) / N_spins;
  specific_heat = energy_variance / (temperature*temperature);
  mean_absolute_magnetization = normalized_expectation_values(4) / N_spins;


  ofile << setprecision(8) << current_cycle;
	ofile << setw(15) << setprecision(8) << temperature;
	ofile << setw(15) << setprecision(8) << mean_energy;
	ofile << setw(15) << setprecision(8) << mean_magnetization;
	ofile << setw(15) << setprecision(8) << specific_heat;
	ofile << setw(15) << setprecision(8) << susceptibility;
	ofile << setw(15) << setprecision(8) << mean_absolute_magnetization;
	ofile << setw(15) << setprecision(8) << accepted_states << "\n";
  ofile.close();
}
