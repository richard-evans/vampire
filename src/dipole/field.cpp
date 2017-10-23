//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>


// Vampire headers
#include "dipole.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "vio.hpp"
#include "errors.hpp"

// dipole module headers
#include "internal.hpp"
#include "material.hpp"
#include "sim.hpp"

namespace dipole{

   //-----------------------------------------------------------------------------
   // Function for updating local temperature fields
   //-----------------------------------------------------------------------------


	void calculate_field(){

      // return if dipole field not enabled
      if(!dipole::activated) return;

		// prevent double calculation for split integration (MPI)
		if(dipole::update_time!=sim::time){

			// Check if update required
		   if(sim::time%dipole::update_rate==0){

			   //if updated record last time at update
			   dipole::update_time=sim::time;

			   // update cell magnetisations
			   cells::mag();
            //MPI::COMM_WORLD.Barrier();
            //fprintf(stderr,"\n >>>> PROBLEMS!!!!!! just after cells::mag()<<<< \n");

			   // recalculate dipole fields
            dipole::internal::update_field();
            //MPI::COMM_WORLD.Barrier();
            //fprintf(stderr,"\n **** PROBLEMS!!!!!! just after dipole::internal::update_field()<<<< \n");

			   // For MPI version, only add local atoms
			   #ifdef MPICF
				   const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
			   #else
				   const int num_local_atoms = dipole::internal::num_atoms;
			   #endif
            //MPI::COMM_WORLD.Barrier();
            //fprintf(stderr,"\n num_local_atoms = %d on my_rank = %d\n",num_local_atoms,vmpi::my_rank);

			   // Update Atomistic Dipolar Field Array
			   for(int atom=0;atom<num_local_atoms;atom++){
				   const int cell = dipole::internal::atom_cell_id_array[atom];
               int type = dipole::internal::atom_type_array[atom];
               //fprintf(stderr,"\tcell = %d x = %f y = %f z = %f mus = %e on my_rank = %d\n",cell,cells::pos_and_mom_array[4*cell+0],cells::pos_and_mom_array[4*cell+1],cells::pos_and_mom_array[4*cell+2],cells::pos_and_mom_array[4*cell+3],vmpi::my_rank);
   	         if(dipole::internal::cells_num_atoms_in_cell[cell]>0 && mp::material[type].non_magnetic==0){
				      // Copy field from macrocell to atomistic spin
				      dipole::atom_dipolar_field_array_x[atom]=dipole::cells_field_array_x[cell];
				      dipole::atom_dipolar_field_array_y[atom]=dipole::cells_field_array_y[cell];
				      dipole::atom_dipolar_field_array_z[atom]=dipole::cells_field_array_z[cell];
                  //fprintf(stderr,"\t atom = %d\tatom_field_x = %f\tatom_field_y = %f\tatom_field_z = %f\ton rank %d\n",atom,dipole::atom_dipolar_field_array_x[atom],dipole::atom_dipolar_field_array_y[atom],dipole::atom_dipolar_field_array_z[atom],vmpi::my_rank);
   	         }
			   }
		   } // End of check for update rate
		} // end of check for update time
   }

} // end of dipole namespace
