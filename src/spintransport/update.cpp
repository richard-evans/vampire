//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spintransport.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

//---------------------------------------------------------------------------------------------------------
// Function to update spin torque field
//---------------------------------------------------------------------------------------------------------
void update(const unsigned int num_local_atoms,            // number of local atoms
            const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
            const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
            const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
            const std::vector<double>& atoms_m_spin_array, // moment of atoms
            std::vector<double>& atoms_x_field_array,      // x-field of atoms
            std::vector<double>& atoms_y_field_array,      // y-field of atoms
            std::vector<double>& atoms_z_field_array       // z-field of atoms
   ){

   //-------------------------------------------------------------------------
   // check that module is needed - if not do nothing
   //-------------------------------------------------------------------------
   if( st::internal::enabled == false ) return;

   //-------------------------------------------------------------------------
   // check that it is time to update
   //-------------------------------------------------------------------------
   if(st::internal::time_counter != st::internal::update_rate){
      // if not increment counter and do nothing
      st::internal::time_counter++;
      return;
   }
   // otherwise reset counter and continue
   else{ st::internal::time_counter = 1; }

   //---------------------------------------------------------------------------------------------------------
   // update cell magnetizations
   //---------------------------------------------------------------------------------------------------------
   st::internal::calculate_cell_magnetization(num_local_atoms, atoms_x_spin_array, atoms_y_spin_array,
                                                               atoms_z_spin_array, atoms_m_spin_array);

   //---------------------------------------------------------------------------------------------------------
   // calculate magnetoresistance
   //---------------------------------------------------------------------------------------------------------
   st::internal::calculate_magnetoresistance();

   //---------------------------------------------------------------------------------------------------------
   // calculate STT field
   //---------------------------------------------------------------------------------------------------------

   return;

}

} // end of spin_transport namespace
