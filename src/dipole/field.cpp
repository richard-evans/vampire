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
#include "vutil.hpp"

// dipole module headers
#include "internal.hpp"
#include "material.hpp"

namespace dipole{

   //-----------------------------------------------------------------------------
   // Function for updating atomic B-field and Hd-field
   //-----------------------------------------------------------------------------


	void calculate_field(const uint64_t sim_time){

      // return if dipole field not enabled
      if(!dipole::activated) return;

		// prevent double calculation for split integration (MPI)
		if(dipole::internal::update_time != sim_time){

			// Check if update required
		   if(sim_time%dipole::update_rate==0){

			   //if updated record last time at update
			   dipole::internal::update_time = sim_time;

            // instantiate timer of cells::mag() function
            //vutil::vtimer_t timer;
            // start timer
            //timer.start();
			   // update cell magnetisations
			   cells::mag();
            // end timer
            //timer.stop();
            // return bandwidth
            //double update_time = timer.elapsed_time();

            //zlog << zTs() << "Calculation cells magnetisation complete. Time taken: " << update_time << "s."<< std::endl;

			   // recalculate dipole fields
            if (dipole::fft == false)  dipole::internal::update_field();
				#ifdef FFT
               else dipole::internal::update_field_fft();
				#endif
					  //MPI::COMM_WORLD.Barrier();
            //fprintf(stderr,"\n **** PROBLEMS!!!!!! just after dipole::internal::update_field()<<<< \n");


			   // For MPI version, only add local atoms
			   #ifdef MPICF
				   const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
			   #else
				   const int num_local_atoms = dipole::internal::num_atoms;
			   #endif

			   // Update Atomistic Dipolar Field and Demag Field Array
			   for(int atom=0;atom<num_local_atoms;atom++){
				   const int cell = dipole::internal::atom_cell_id_array[atom];
               int type = dipole::internal::atom_type_array[atom];
   	         if(dipole::internal::cells_num_atoms_in_cell[cell]>0 && mp::material[type].non_magnetic==0){
				      // Copy B-field from macrocell to atomistic spin
				      dipole::atom_dipolar_field_array_x[atom]=dipole::cells_field_array_x[cell];
				      dipole::atom_dipolar_field_array_y[atom]=dipole::cells_field_array_y[cell];
				      dipole::atom_dipolar_field_array_z[atom]=dipole::cells_field_array_z[cell];
                  // Unroll Hdemag field
				      dipole::atom_mu0demag_field_array_x[atom]=dipole::cells_mu0Hd_field_array_x[cell];
				      dipole::atom_mu0demag_field_array_y[atom]=dipole::cells_mu0Hd_field_array_y[cell];
				      dipole::atom_mu0demag_field_array_z[atom]=dipole::cells_mu0Hd_field_array_z[cell];
   	         }
			   }
		   } // End of check for update rate
		} // end of check for update time
   }

} // end of dipole namespace
