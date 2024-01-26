// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdio.h>

// Vampire Header files
#include "LSF_RK4.hpp"
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "internal.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Field calculation functions
int calculate_spin_fields(const int, const int);
int calculate_external_fields(const int, const int);

namespace LSF_RK4_arrays{

   // Local arrays for LSF integration
   std::vector<double> x_lsf_array;
   std::vector<double> y_lsf_array;
   std::vector<double> z_lsf_array;

   std::vector<double> x_initial_spin_array;
   std::vector<double> y_initial_spin_array;
   std::vector<double> z_initial_spin_array;

   std::vector<double> x_k1_array;
	std::vector<double> y_k1_array;
	std::vector<double> z_k1_array;

	std::vector<double> x_k2_array;
	std::vector<double> y_k2_array;
	std::vector<double> z_k2_array;

	std::vector<double> x_k3_array;
	std::vector<double> y_k3_array;
	std::vector<double> z_k3_array;

	std::vector<double> x_k4_array;
	std::vector<double> y_k4_array;
	std::vector<double> z_k4_array;

   // Flag to define state of LSF arrays (initialised/uninitialised)
   bool LSF_RK4_set = false;

   std::vector<double> mod_S;

   std::vector<double> tx;
   std::vector<double> ty;
   std::vector<double> tz;

}
namespace sim{

   int LSF_RK4_init(){
      // Check calling of routine if error checking is activated
      if (err::check == true){
         std::cout << "sim:LSF_RK4_init has been called" << std::endl;
      }

      using namespace LSF_RK4_arrays;

      x_lsf_array.resize(atoms::num_atoms, 0.0);
      y_lsf_array.resize(atoms::num_atoms, 0.0);
      z_lsf_array.resize(atoms::num_atoms, 0.0);

      x_initial_spin_array.resize(atoms::num_atoms, 0.0);
      y_initial_spin_array.resize(atoms::num_atoms, 0.0);
      z_initial_spin_array.resize(atoms::num_atoms, 0.0);

		x_k1_array.resize(atoms::num_atoms, 0.0);
		y_k1_array.resize(atoms::num_atoms, 0.0);
		z_k1_array.resize(atoms::num_atoms, 0.0);

		x_k2_array.resize(atoms::num_atoms, 0.0);
		y_k2_array.resize(atoms::num_atoms, 0.0);
		z_k2_array.resize(atoms::num_atoms, 0.0);

		x_k3_array.resize(atoms::num_atoms, 0.0);
		y_k3_array.resize(atoms::num_atoms, 0.0);
		z_k3_array.resize(atoms::num_atoms, 0.0);

		x_k4_array.resize(atoms::num_atoms, 0.0);
		y_k4_array.resize(atoms::num_atoms, 0.0);
		z_k4_array.resize(atoms::num_atoms, 0.0);

      mod_S.resize(atoms::num_atoms, 1.0);

      tx.resize(atoms::num_atoms, 0.0);
      ty.resize(atoms::num_atoms, 0.0);
      tz.resize(atoms::num_atoms, 0.0);

      // Disable external thermal field calculations
      sim::hamiltonian_simulation_flags[3] = 0;

      LSF_RK4_set = true;

      return EXIT_SUCCESS;
   }

   // LSF-RK4 magnetic field function
   void calculate_lsf_rk4_magnetic_field(const int start_index, const int end_index){

      // Check calling of routine if error checking is activated
      if (err::check == true){
         std::cout << "calculate_lsf_rk4_magnetic_field has been called" << std::endl;
      }

      // LSF-RK4 Hamiltonian calculation
      for (int atom = start_index; atom < end_index; atom++){

         const int imaterial = atoms::type_array[atom];

         const double sx = atoms::x_spin_array[atom];
         const double sy = atoms::y_spin_array[atom];
         const double sz = atoms::z_spin_array[atom];

         const double imu_S = -1.0 / mp::material[imaterial].mu_s_SI;

         const double L2 = 2.0 * sim::internal::lsf_second_order_coefficient[imaterial] * imu_S;
         const double L4 = 4.0 * sim::internal::lsf_fourth_order_coefficient[imaterial] * imu_S;
         const double L6 = 6.0 * sim::internal::lsf_sixth_order_coefficient[imaterial] * imu_S;

         const double ss2 = sx*sx + sy*sy + sz*sz;

         LSF_RK4_arrays::x_lsf_array[atom] = L2*sx + L4*sx*ss2 + L6*sx*ss2*ss2;
         LSF_RK4_arrays::y_lsf_array[atom] = L2*sy + L4*sy*ss2 + L6*sy*ss2*ss2;
         LSF_RK4_arrays::z_lsf_array[atom] = L2*sz + L4*sz*ss2 + L6*sz*ss2*ss2;
      }
   }

   namespace internal{

      void lsf_rk4_step(){

         // Check calling of routine if error checking is activated
         if (err::check == true){
            std::cout << "sim::LSF_RK4 has been called" << std::endl;
         }

         using namespace LSF_RK4_arrays;

         // Check for initialisation of LSF integration arrays
         if (LSF_RK4_set == false)
            sim::LSF_RK4_init();

         // Local variables for system integration
         const int num_atoms = atoms::num_atoms;
         double xyz[3];   // Local delta spin components
         double S_new[3]; // New local spin moment
         const double kB = 1.3806503e-23;

         // Calculate fields
         calculate_spin_fields(0, num_atoms);
         calculate_lsf_rk4_magnetic_field(0, num_atoms);
         calculate_external_fields(0, num_atoms);

         // Store initial spin positions
         for (int atom = 0; atom < num_atoms; atom++){
            x_initial_spin_array[atom] = atoms::x_spin_array[atom];
            y_initial_spin_array[atom] = atoms::y_spin_array[atom];
            z_initial_spin_array[atom] = atoms::z_spin_array[atom];
         }

         // Thermal noise based on Gaussian function
         double sigma = (sqrt((2.0 * kB * sim::temperature * mp::gamma_SI) / (mp::dt_SI)));
         generate(tx.begin(), tx.begin() + num_atoms, mtrandom::gaussian);
         generate(ty.begin(), ty.begin() + num_atoms, mtrandom::gaussian);
         generate(tz.begin(), tz.begin() + num_atoms, mtrandom::gaussian);

         // Calculate K1
         for (int atom = 0; atom < num_atoms; atom++){

            const int imaterial = atoms::type_array[atom];
            const double alpha = mp::material[atoms::type_array[atom]].alpha;
            const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

            // Store local spin in S and local field in H
            const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};

            const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LSF_RK4_arrays::x_lsf_array[atom],
                                 atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LSF_RK4_arrays::y_lsf_array[atom],
                                 atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LSF_RK4_arrays::z_lsf_array[atom]};

            // Calculate Delta S
            xyz[0] = (-mp::gamma_SI * (S[1] * H[2] - S[2] * H[1])) + (mp::gamma_SI * alpha * H[0]) + (tx[atom] * (sigma * sqrt(alpha / mu)));
            xyz[1] = (-mp::gamma_SI * (S[2] * H[0] - S[0] * H[2])) + (mp::gamma_SI * alpha * H[1]) + (ty[atom] * (sigma * sqrt(alpha / mu)));
            xyz[2] = (-mp::gamma_SI * (S[0] * H[1] - S[1] * H[0])) + (mp::gamma_SI * alpha * H[2]) + (tz[atom] * (sigma * sqrt(alpha / mu)));

				// k1 = \Delta S * \Delta t
            x_k1_array[atom] = xyz[0] * mp::dt_SI;
            y_k1_array[atom] = xyz[1] * mp::dt_SI;
            z_k1_array[atom] = xyz[2] * mp::dt_SI;

            // Write new spin vectors
            atoms::x_spin_array[atom] += 0.5*x_k1_array[atom];
            atoms::y_spin_array[atom] += 0.5*y_k1_array[atom];
				atoms::z_spin_array[atom] += 0.5*z_k1_array[atom];

         }

         // Recalculate spin dependent fields
         calculate_spin_fields(0, num_atoms);
         calculate_lsf_rk4_magnetic_field(0, num_atoms);

         // Calculate K2
         for (int atom = 0; atom < num_atoms; atom++){

            const int imaterial = atoms::type_array[atom];
            const double alpha = mp::material[atoms::type_array[atom]].alpha;
            const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

            // Store local spin in S and local field in H
            const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};

            const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LSF_RK4_arrays::x_lsf_array[atom],
                                 atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LSF_RK4_arrays::y_lsf_array[atom],
                                 atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LSF_RK4_arrays::z_lsf_array[atom]};

            // Calculate Delta S
            xyz[0] = (-mp::gamma_SI * (S[1] * H[2] - S[2] * H[1])) + (mp::gamma_SI * alpha * H[0]) + (tx[atom] * (sigma * sqrt(alpha / mu)));
            xyz[1] = (-mp::gamma_SI * (S[2] * H[0] - S[0] * H[2])) + (mp::gamma_SI * alpha * H[1]) + (ty[atom] * (sigma * sqrt(alpha / mu)));
            xyz[2] = (-mp::gamma_SI * (S[0] * H[1] - S[1] * H[0])) + (mp::gamma_SI * alpha * H[2]) + (tz[atom] * (sigma * sqrt(alpha / mu)));

            // k2
            x_k2_array[atom] = xyz[0] * mp::dt_SI;
            y_k2_array[atom] = xyz[1] * mp::dt_SI;
            z_k2_array[atom] = xyz[2] * mp::dt_SI;

            // Write new spin vectors
            atoms::x_spin_array[atom] += 0.5*x_k2_array[atom];
            atoms::y_spin_array[atom] += 0.5*y_k2_array[atom];
				atoms::z_spin_array[atom] += 0.5*z_k2_array[atom];
         }

			// Calculate K3
         for (int atom = 0; atom < num_atoms; atom++){

            const int imaterial = atoms::type_array[atom];
            const double alpha = mp::material[atoms::type_array[atom]].alpha;
            const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

            // Store local spin in S and local field in H
            const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};

            const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LSF_RK4_arrays::x_lsf_array[atom],
                                 atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LSF_RK4_arrays::y_lsf_array[atom],
                                 atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LSF_RK4_arrays::z_lsf_array[atom]};

            // Calculate Delta S
            xyz[0] = (-mp::gamma_SI * (S[1] * H[2] - S[2] * H[1])) + (mp::gamma_SI * alpha * H[0]) + (tx[atom] * (sigma * sqrt(alpha / mu)));
            xyz[1] = (-mp::gamma_SI * (S[2] * H[0] - S[0] * H[2])) + (mp::gamma_SI * alpha * H[1]) + (ty[atom] * (sigma * sqrt(alpha / mu)));
            xyz[2] = (-mp::gamma_SI * (S[0] * H[1] - S[1] * H[0])) + (mp::gamma_SI * alpha * H[2]) + (tz[atom] * (sigma * sqrt(alpha / mu)));

            // k3
            x_k3_array[atom] = xyz[0] * mp::dt_SI;
            y_k3_array[atom] = xyz[1] * mp::dt_SI;
            z_k3_array[atom] = xyz[2] * mp::dt_SI;

            // Write new spin vectors
            atoms::x_spin_array[atom] += x_k3_array[atom];
            atoms::y_spin_array[atom] += y_k3_array[atom];
				atoms::z_spin_array[atom] += z_k3_array[atom];

         }

			// Recalculate spin dependent fields
         calculate_spin_fields(0, num_atoms);
         calculate_lsf_rk4_magnetic_field(0, num_atoms);

			// Calculate K4
         for (int atom = 0; atom < num_atoms; atom++){

            const int imaterial = atoms::type_array[atom];
            const double alpha = mp::material[atoms::type_array[atom]].alpha;
            const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

            // Store local spin in S and local field in H
            const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};

            const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LSF_RK4_arrays::x_lsf_array[atom],
                                 atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LSF_RK4_arrays::y_lsf_array[atom],
                                 atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LSF_RK4_arrays::z_lsf_array[atom]};

            // Calculate Delta S
            xyz[0] = (-mp::gamma_SI * (S[1] * H[2] - S[2] * H[1])) + (mp::gamma_SI * alpha * H[0]) + (tx[atom] * (sigma * sqrt(alpha / mu)));
            xyz[1] = (-mp::gamma_SI * (S[2] * H[0] - S[0] * H[2])) + (mp::gamma_SI * alpha * H[1]) + (ty[atom] * (sigma * sqrt(alpha / mu)));
            xyz[2] = (-mp::gamma_SI * (S[0] * H[1] - S[1] * H[0])) + (mp::gamma_SI * alpha * H[2]) + (tz[atom] * (sigma * sqrt(alpha / mu)));

            // k4
            x_k4_array[atom] = xyz[0] * mp::dt_SI;
            y_k4_array[atom] = xyz[1] * mp::dt_SI;
            z_k4_array[atom] = xyz[2] * mp::dt_SI;

         }

         // Calculate RK4 Step
         for (int atom = 0; atom < num_atoms; atom++){
            atoms::x_spin_array[atom] = x_initial_spin_array[atom] + (1.0/6.0)*(x_k1_array[atom] + (2.0*x_k2_array[atom]) + (2.0*x_k3_array[atom]) + x_k4_array[atom]);
            atoms::y_spin_array[atom] = y_initial_spin_array[atom] + (1.0/6.0)*(y_k1_array[atom] + (2.0*y_k2_array[atom]) + (2.0*y_k3_array[atom]) + y_k4_array[atom]);
            atoms::z_spin_array[atom] = z_initial_spin_array[atom] + (1.0/6.0)*(z_k1_array[atom] + (2.0*z_k2_array[atom]) + (2.0*z_k3_array[atom]) + z_k4_array[atom]);
         }

         // Store spin length data
         for (int atom = 0; atom < num_atoms; atom++){
            const double sx = atoms::x_spin_array[atom];
            const double sy = atoms::y_spin_array[atom];
            const double sz = atoms::z_spin_array[atom];

            mod_S[atom] = sqrt(sx * sx + sy * sy + sz * sz);
         }

         return;
      }

   } // End of indental namespace

} // End of sim namespace
