//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>
#include "vio.hpp"

// Vampire headers
#include "spinpumping.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_pumping{

//---------------------------------------------------------------------------------------------------------
// Function to update spin torque field
//---------------------------------------------------------------------------------------------------------
void update(const unsigned int num_local_atoms,            // number of local atoms
            const uint64_t time_sim,                       // simulation time 
            const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
            const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
            const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
            const std::vector<double>& atoms_m_spin_array  // moment of atoms
   ){

   //-------------------------------------------------------------------------
   // check that module is needed - if not do nothing
   //-------------------------------------------------------------------------
   if( spin_pumping::internal::enabled == false ) return;

   //-------------------------------------------------------------------------
   // check that it is time to update
   //-------------------------------------------------------------------------
   if(spin_pumping::internal::time_counter != spin_pumping::internal::update_rate){
      // if not increment counter and do nothing
      spin_pumping::internal::time_counter++;
      return;
   }
   // otherwise reset counter and continue
   else{ 
      spin_pumping::internal::time_counter = 1; 
      spin_pumping::internal::config_counter += 1; 
   }

   //---------------------------------------------------------------------------------------------------------
   // calculate atomistic spin pumping
   //---------------------------------------------------------------------------------------------------------
   // If enabled, output calculated atomistic coordinates and moments (passing local values)

   // Get previous spin configuration
   const std::vector <double> atoms_x_old_spin_array = spin_pumping::internal::get_old_spins_x(num_local_atoms);
   const std::vector <double> atoms_y_old_spin_array = spin_pumping::internal::get_old_spins_y(num_local_atoms);
   const std::vector <double> atoms_z_old_spin_array = spin_pumping::internal::get_old_spins_z(num_local_atoms);
   spin_pumping::internal::calculate_spin_pumping(num_local_atoms, time_sim, atoms_x_spin_array, atoms_y_spin_array,
                                                   atoms_z_spin_array, atoms_x_old_spin_array, atoms_y_old_spin_array,
                                                   atoms_z_old_spin_array, atoms_m_spin_array);

   //---------------------------------------------------------------------------------------------------------
   // calculate cells spin pumping
   //---------------------------------------------------------------------------------------------------------
   // spin_pumping::internal::calculate_cell_magnetization(num_local_atoms, atoms_x_spin_array, atoms_y_spin_array,
   //                                                             atoms_z_spin_array, atoms_m_spin_array);
   // spin_pumping::internal::calculate_cells_spin_pumping();

   // update counter for spin pumping output
   const uint64_t tmp_counter = spin_pumping::internal::config_counter-1;
   // Output atoms only
   if(spin_pumping::internal::output_atomistic_spin_pumping_flag){ spin_pumping::internal::output_atomistic_spin_pumping(tmp_counter);}
   // Output cells only
   if(spin_pumping::internal::output_cells_spin_pumping_flag){ spin_pumping::internal::output_atomistic_spin_pumping(tmp_counter);}

   return;

}

} // end of spin_pumping namespace
