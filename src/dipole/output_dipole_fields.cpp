//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2025. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
//#include <cmath>
#include <iostream>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
//#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// dipole module headers
#include "internal.hpp"

// alias interal dipole namespace for brevity
namespace dp = dipole::internal;

namespace dipole{

namespace internal{

//------------------------------------------------------------------------------
// Function to output calculated dipole cell coords and fields on root process
//------------------------------------------------------------------------------
void output_dipole_fields(){

   // If dipole output not enabled, then do nothing
   if(!dipole::internal::output_dipole_field) return;

   // inform user that dipole fields are being outputted
   zlog << zTs() << "Outputting dipole fields to file" << std::endl;

   // if rank = 0 open output file
   if(vmpi::my_rank == 0){

      dp_fields.open("dipole-field.txt");

      for (int i = 0 ; i < dipole::internal::cells_num_cells; i ++){
         dp_fields << i << "\t" <<  // cell ID
         cells::pos_and_mom_array[4*i+0] << "\t" << // x
         cells::pos_and_mom_array[4*i+1] << "\t" << // y
         cells::pos_and_mom_array[4*i+2] << "\t" << // z
         cells::pos_and_mom_array[4*i+3] << "\t" << // m
         dipole::cells_field_array_x[i] << "\t" << // Bx
         dipole::cells_field_array_y[i] << "\t" << // By
         dipole::cells_field_array_z[i] << "\t" << std::endl; // Bz
      }

   }

   return;

}

} // end of namespace internal
} // end of namespace dipole
