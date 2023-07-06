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
#include "gpu.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "vio.hpp"
#include "errors.hpp"
#include "vutil.hpp"

// dipole module headers
#include "internal.hpp"
#include "material.hpp"
#include "../hierarchical/internal.hpp"
#include "hierarchical.hpp"


namespace dipole{

   namespace internal{
      void calculate_macrocell_dipole_field();
   }

   //-----------------------------------------------------------------------------
   // Function for updating atomic B-field and Hd-field
   //-----------------------------------------------------------------------------
	void calculate_field(const uint64_t sim_time,
                        std::vector <double>& x_spin_array, // atomic spin directions
                        std::vector <double>& y_spin_array,
                        std::vector <double>& z_spin_array,
                        std::vector <double>& m_spin_array, // atomic spin moment
                        std::vector < bool >& magnetic){ // is magnetic

      // return if dipole field not enabled
      if(!dipole::activated) return;

		// prevent double calculation for split integration (MPI)
		if(dipole::internal::update_time != static_cast<int>(sim_time)){

			// Check if update required
		   if(sim_time%dipole::update_rate == 0){

			   //if updated record last time at update
			   dipole::internal::update_time = sim_time;

            // // for gpu acceleration, transfer spin positions now (does nothing for serial)
            // gpu::transfer_spin_positions_from_gpu_to_cpu();

            switch (dipole::internal::solver){

               case dipole::internal::macrocell:
                  dipole::internal::calculate_macrocell_dipole_field();
                  break;

               case dipole::internal::tensor:
                  #ifdef CUDA
               	   gpu::update_dipolar_fields();
                  #else
                     dipole::internal::calculate_macrocell_dipole_field();
                  #endif
                  break;

               case dipole::internal::atomistic:
                  dipole::internal::calculate_atomistic_dipole_field(x_spin_array, y_spin_array, z_spin_array);
                  break;

               case dipole::internal::hierarchical:
                  hierarchical::update(x_spin_array, y_spin_array, z_spin_array, m_spin_array, magnetic);
                  break;

               case dipole::internal::fft:
                  dipole::internal::update_field_fft();
                  break;

               case dipole::internal::atomisticfft:
                  dipole::internal::atomistic_fft::update_field_atomistic_fft();
                  break;


            }

            // // for gpu acceleration, transfer calculated fields now (does nothing for serial)
            // gpu::transfer_dipole_fields_from_cpu_to_gpu();
            // // for gpu acceleration, transfer calculated cells dipolar fields now (does nothing for serial)
            // gpu::transfer_dipole_cells_fields_from_gpu_to_cpu();

		   } // End of check for update rate
		} // end of check for update time

      return;

   }

   namespace internal{

      void calculate_macrocell_dipole_field(){
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
         dipole::internal::update_field();

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

            if(dipole::internal::cells_num_atoms_in_cell[cell]>0 && mp::material[type].non_magnetic==false){

               // Copy B-field from macrocell to atomistic spin
               // Copy B-field from macrocell to atomistic spin
               dipole::atom_dipolar_field_array_x[atom] = dipole::cells_field_array_x[cell];
               dipole::atom_dipolar_field_array_y[atom] = dipole::cells_field_array_y[cell];
               dipole::atom_dipolar_field_array_z[atom] = dipole::cells_field_array_z[cell];

               // Unroll Hdemag field
               dipole::atom_mu0demag_field_array_x[atom] = dipole::cells_mu0Hd_field_array_x[cell];
               dipole::atom_mu0demag_field_array_y[atom] = dipole::cells_mu0Hd_field_array_y[cell];
               dipole::atom_mu0demag_field_array_z[atom] = dipole::cells_mu0Hd_field_array_z[cell];

            }
         }

         return;

      } // end of function

} // end of internal namespace

} // end of dipole namespace
