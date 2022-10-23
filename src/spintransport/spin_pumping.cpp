//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo & Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "spintransport.hpp"
#include "material.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{
	namespace internal{

   	//---------------------------------------------------------------------------
   	// Function to compute atomistic s_i x ds_i/dt
   	//---------------------------------------------------------------------------
		void calculate_spin_cross_spin_time_derivative(const unsigned int num_local_atoms,            // number of local atoms
		                                  const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
		                                  const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
		                                  const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
		                                  const std::vector<double>& atoms_x_old_spin_array, // old x-spin vector of atoms
		                                  const std::vector<double>& atoms_y_old_spin_array, // old y-spin vector of atoms
		                                  const std::vector<double>& atoms_z_old_spin_array, // old z-spin-vector of atoms
		                                  const std::vector<double>& atoms_m_spin_array  // moment of atoms
		                               ){
				// Arrays to store calculation ds/dt
				std::vector <double> x_dsdt_array(num_local_atoms,0.0);
				std::vector <double> y_dsdt_array(num_local_atoms,0.0);
				std::vector <double> z_dsdt_array(num_local_atoms,0.0);

				// Calculate time derivative of spin
   			for(unsigned int atom = 0; atom < num_local_atoms; atom++){
					const double dt_inv = 1.0/mp::dt;
      			x_dsdt_array[atom] = (atoms_x_spin_array[atom]-atoms_x_old_spin_array[atom])*dt_inv;
      			y_dsdt_array[atom] = (atoms_y_spin_array[atom]-atoms_y_old_spin_array[atom])*dt_inv;
      			z_dsdt_array[atom] = (atoms_z_spin_array[atom]-atoms_z_old_spin_array[atom])*dt_inv;
   			}

				// Calculate cross-product derivative of spin
   			for(unsigned int atom = 0; atom < num_local_atoms; atom++){
   			   // get magnetic moment (muB)
   			   const double mm_sq = atoms_m_spin_array[atom]*atoms_m_spin_array[atom];
   			   const double mm_sq_inv = 1.0/atoms_m_spin_array[atom]*atoms_m_spin_array[atom];
					const double  a_x = atoms_x_spin_array[atom];
					const double  a_y = atoms_y_spin_array[atom];
					const double  a_z = atoms_z_spin_array[atom];
					const double  b_x = x_dsdt_array[atom];
					const double  b_y = y_dsdt_array[atom];
					const double  b_z = z_dsdt_array[atom];
         		st::internal::x_s_cross_dsdt_array[atom] = (a_y * b_z - a_z * b_y)/**mm_sq_inv*/;
         		st::internal::y_s_cross_dsdt_array[atom] = (a_z * b_x - a_x * b_z)/**mm_sq_inv*/;
         		st::internal::z_s_cross_dsdt_array[atom] = (a_x * b_y - a_y * b_x)/**mm_sq_inv*/;
   			}
													
		} // end of function
	} // end of namespace
} // end of namespace
