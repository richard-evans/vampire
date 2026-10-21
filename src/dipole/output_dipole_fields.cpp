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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "material.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// dipole module headers
#include "internal.hpp"

// alias interal dipole namespace for brevity
namespace dp = dipole::internal;

namespace dipole{

//------------------------------------------------------------------------------
// Function to output calculated dipole cell coords and fields on root process
// (externally visible: called from both the tensor/macrocell and hierarchical
// dipole solvers, so it lives in the public dipole:: namespace rather than
// dipole::internal::)
//------------------------------------------------------------------------------
void output_dipole_fields(){

   // If dipole output not enabled, then do nothing
   if(!dipole::internal::output_dipole_field) return;

   // inform user that dipole fields are being outputted
   zlog << zTs() << "Outputting dipole fields to file" << std::endl;

   // Calculate total num atoms in cell (across the whole system)
   std::vector<int> num_atoms_out = cells::num_atoms_in_cell;
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_out[0], dipole::cells_num_cells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // The atomistic solver does not calculate cell fields, so use the mean field of the magnetic atoms in each cell
   const bool atomistic_solver = dipole::internal::solver == dipole::internal::atomistic;
   std::vector<double> mean_field_x(0);
   std::vector<double> mean_field_y(0);
   std::vector<double> mean_field_z(0);
   if(atomistic_solver){
      mean_field_x.assign(dipole::cells_num_cells, 0.0);
      mean_field_y.assign(dipole::cells_num_cells, 0.0);
      mean_field_z.assign(dipole::cells_num_cells, 0.0);
      std::vector<double> num_magnetic_atoms(dipole::cells_num_cells, 0.0);
      const int num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      for(int atom = 0; atom < num_local_atoms; atom++){
         if(mp::material[dp::atom_type_array[atom]].non_magnetic) continue;
         const int cell = dp::atom_cell_id_array[atom];
         mean_field_x[cell] += dipole::atom_dipolar_field_array_x[atom];
         mean_field_y[cell] += dipole::atom_dipolar_field_array_y[atom];
         mean_field_z[cell] += dipole::atom_dipolar_field_array_z[atom];
         num_magnetic_atoms[cell] += 1.0;
      }
      vmpi::all_reduce_sum(mean_field_x);
      vmpi::all_reduce_sum(mean_field_y);
      vmpi::all_reduce_sum(mean_field_z);
      vmpi::all_reduce_sum(num_magnetic_atoms);
      for(int i = 0; i < dipole::cells_num_cells; i++){
         if(num_magnetic_atoms[i] > 0.0){
            mean_field_x[i] /= num_magnetic_atoms[i];
            mean_field_y[i] /= num_magnetic_atoms[i];
            mean_field_z[i] /= num_magnetic_atoms[i];
         }
      }
   }
   const std::vector<double>& field_x = atomistic_solver ? mean_field_x : dipole::cells_field_array_x;
   const std::vector<double>& field_y = atomistic_solver ? mean_field_y : dipole::cells_field_array_y;
   const std::vector<double>& field_z = atomistic_solver ? mean_field_z : dipole::cells_field_array_z;

   // if rank = 0 open output file
   if(vmpi::my_rank == 0){

      std::ofstream dp_fields("dipole-field.txt");
      dp_fields << std::setprecision(16);

      for (int i = 0 ; i < dipole::cells_num_cells; i ++){
         dp_fields << i << "\t" <<  // cell ID
         num_atoms_out[i] << "\t" << // num atoms in cell (globally summed, see comment above)
         cells::pos_and_mom_array[4*i+0] << "\t" << // x
         cells::pos_and_mom_array[4*i+1] << "\t" << // y
         cells::pos_and_mom_array[4*i+2] << "\t" << // z
         cells::pos_and_mom_array[4*i+3] << "\t" << // m
         field_x[i] << "\t" << // Bx
         field_y[i] << "\t" << // By
         field_z[i] << "\t" << std::endl; // Bz
      }

   }

   return;

}

} // end of namespace dipole
