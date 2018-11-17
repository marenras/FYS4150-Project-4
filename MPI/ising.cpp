#include "ising.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <random>

Ising::Ising(int dim_lattice)
{
  // Setting "self" values
  this->dim_lattice = dim_lattice;

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
}


void Ising::InitializeLattice(double temperature, bool random_lattice)
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

      if (random_lattice == false)
      {
        lattice(x,y) = 1.0;
      }
      else
      {
          double rand_condition = RNG(generator);
          if (rand_condition < 0.5) {lattice(x,y) = 1.0;}
          else {lattice(x,y) = -1.0;}
      }
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
         lattice(rand_x, PBC(rand_y, dim_lattice, 1))  +
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
