//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo & Richard F L Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "constants.hpp"
#include "material.hpp"
#include "spinpumping.hpp"

// spinpumping module headers
#include "internal.hpp"

namespace spin_pumping{
	namespace internal{

   	//---------------------------------------------------------------------------
   	// Function to compute atomistic spin pumping as: s_i x ds_i/dt
   	//---------------------------------------------------------------------------
		void calculate_spin_pumping(const unsigned int num_local_atoms,            // number of local atoms
		                                  const uint64_t time_sim, 						// simulation time 
		                                  const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
		                                  const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
		                                  const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
		                                  const std::vector<double>& atoms_x_old_spin_array, // old x-spin vector of atoms
		                                  const std::vector<double>& atoms_y_old_spin_array, // old y-spin vector of atoms
		                                  const std::vector<double>& atoms_z_old_spin_array, // old z-spin-vector of atoms
		                                  const std::vector<double>& atoms_m_spin_array  // moment of atoms
		                               ){

				//-------------------------------------------------------------------------
				// check that is not first step - if not do nothing
				//-------------------------------------------------------------------------
				if( time_sim-1<1 ){return;}

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

            // Calculate prefactor h/e^2 to get spin current in correct units
            const double h_o_e2 = 8.9674108e-35;  // in J s; to pass from spin-mixing onductance in in 1/(Ohm m^2) to effective spin-mixing conductance in 1/m^2
            const double hbar_o_2muB = 0.5 * (1.05457182e-34/constants::muB);  // to convert spin polarisation into spin current
            const double prefactor = h_o_e2 * hbar_o_2muB;
            // Calculate cross-product derivative of spin
            for(unsigned int atom = 0; atom < num_local_atoms; atom++){
               const int mat = spin_pumping::internal::atoms_type_array[atom];
               if(spin_pumping::internal::material_magnetic[mat]){
                  // get magnetic moment (muB)
                  const double mus = mp::material[mat].mu_s_SI;
                  const double mm_sq = atoms_m_spin_array[atom]*atoms_m_spin_array[atom];
                  const double one_o_mm_sq = 1.0/mm_sq;
                  const double spin_mix_conductance = spin_pumping::internal::mp[mat].spin_mix_conductance.get();
                  const double prefactor_atom = prefactor * spin_mix_conductance * one_o_mm_sq;
                  const double  a_x = atoms_x_spin_array[atom];
                  const double  a_y = atoms_y_spin_array[atom];
                  const double  a_z = atoms_z_spin_array[atom];
                  const double  b_x = x_dsdt_array[atom];
                  const double  b_y = y_dsdt_array[atom];
                  const double  b_z = z_dsdt_array[atom];
                  spin_pumping::internal::x_atom_spin_pumping_array[atom] = /*prefactor_atom **/ (a_y * b_z - a_z * b_y);
                  spin_pumping::internal::y_atom_spin_pumping_array[atom] = /*prefactor_atom **/ (a_z * b_x - a_x * b_z);
                  spin_pumping::internal::z_atom_spin_pumping_array[atom] = /*prefactor_atom **/ (a_x * b_y - a_y * b_x);
   			   }
   	      }

		   } // end of function

   	   //---------------------------------------------------------------------------
   	   // Function to compute atomistic spin pumping as: s_i x ds_i/dt
   	   //---------------------------------------------------------------------------
		   // void calculate_cells_spin_pumping(
         // } // end of function

      } // end of namespace
} // end of namespace
