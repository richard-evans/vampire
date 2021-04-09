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

   //---------------------------------------------------------------------------
   // Function to calculate spin transfer torque field for each atom
   //---------------------------------------------------------------------------
   void calculate_field(const unsigned int start_index,            // first atom
                        const unsigned int end_index,              // last atom
                        std::vector<double>& atoms_x_field_array,  // x-field of atoms
                        std::vector<double>& atoms_y_field_array,  // y-field of atoms
                        std::vector<double>& atoms_z_field_array   // z-field of atoms
      ){

         //-------------------------------------------------------------------------
         // check that module is needed - if not do nothing
         //-------------------------------------------------------------------------
         if( st::internal::enabled == false ) return;

         //---------------------------------------------------------------------------
         // loop over all atoms and apply cell spin torque field
         //---------------------------------------------------------------------------
         for(unsigned int atom = start_index; atom < end_index; atom++){

            // get cell id
            const uint64_t cell = st::internal::atom_in_cell[atom];

            atoms_x_field_array[atom] += st::internal::cell_spin_torque_fields[3*cell+0];
            atoms_y_field_array[atom] += st::internal::cell_spin_torque_fields[3*cell+1];
            atoms_z_field_array[atom] += st::internal::cell_spin_torque_fields[3*cell+2];

         }

      return;

   }

}
