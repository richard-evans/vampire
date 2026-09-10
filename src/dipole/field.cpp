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
      void broadcast_cell_field_to_atoms();

      // Human-readable name for the active solver, for scaling-study log output
      const char* dipole_solver_name(){
         switch (dipole::internal::solver){
            case dipole::internal::macrocell:   return "macrocell";
            case dipole::internal::tensor:      return "tensor";
            case dipole::internal::atomistic:   return "atomistic";
            case dipole::internal::hierarchical:return "hierarchical";
            case dipole::internal::fft:         return "fft";
            case dipole::internal::atomisticfft:return "atomisticfft";
            default:                            return "unknown";
         }
      }
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

            // Time each individual dipole field recalculation, for the
            // ARCHER2 weak/strong scaling study (papers/dipole). Barrier
            // first so slow ranks don't make fast ranks look slow, then
            // report the slowest rank's time (the one that sets wall clock).
            vmpi::barrier();
            vutil::vtimer_t dipole_scaling_timer;
            dipole_scaling_timer.start();

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
                  // hierarchical::update() only fills the per-macrocell
                  // dipole::cells_field_array_*; broadcast it down to
                  // dipole::atom_dipolar_field_array_* (read by
                  // dipole:output-atomistic-dipole-field) the same way the
                  // macrocell/tensor path does below.
                  dipole::internal::broadcast_cell_field_to_atoms();
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

            dipole_scaling_timer.stop();
            const double dipole_scaling_local_time = dipole_scaling_timer.elapsed_time();
            const double dipole_scaling_max_time = vmpi::reduce_max(dipole_scaling_local_time);
            if(vmpi::master){
               std::cout << "DIPOLE_SCALING UPDATE solver=" << dipole::internal::dipole_solver_name()
                          << " nprocs=" << vmpi::num_processors
                          << " sim_time=" << sim_time
                          << " update_time_s=" << dipole_scaling_max_time << std::endl;
               zlog << zTs() << "DIPOLE_SCALING UPDATE solver=" << dipole::internal::dipole_solver_name()
                    << " nprocs=" << vmpi::num_processors
                    << " sim_time=" << sim_time
                    << " update_time_s=" << dipole_scaling_max_time << std::endl;
            }

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

         // Update Atomistic Dipolar Field and Demag Field Array
         broadcast_cell_field_to_atoms();

         return;

      } // end of function

      //------------------------------------------------------------------------
      // Copy each atom's host macrocell dipole (and demag) field down onto
      // dipole::atom_dipolar_field_array_*/atom_mu0demag_field_array_*. Shared
      // by the macrocell/tensor path (calculate_macrocell_dipole_field(), above)
      // and the hierarchical path (calculate_field(), which has no equivalent
      // broadcast of its own since hierarchical::update() only fills the
      // per-cell field arrays).
      //
      // Non-magnetic atoms (e.g. a dense non-magnetic sensor/read-back layer
      // stacked above a granular medium) are deliberately included here: they
      // carry no moment and so contribute nothing to cells::mag() or to any
      // LLG torque, but they still sit in a macrocell and are valid dipole
      // field *probes* for dipole:output-atomistic-dipole-field. Only skip a
      // cell with no atoms in it at all (field undefined/meaningless there).
      //------------------------------------------------------------------------
      void broadcast_cell_field_to_atoms(){

         // For MPI version, only add local atoms
         #ifdef MPICF
            const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
         #else
            const int num_local_atoms = dipole::internal::num_atoms;
         #endif

         for(int atom=0;atom<num_local_atoms;atom++){

            const int cell = dipole::internal::atom_cell_id_array[atom];

            // No guard on dipole::internal::cells_num_atoms_in_cell[cell] here:
            // that count is magnetic atoms only, and would incorrectly skip
            // every atom (magnetic or not) hosted in a non-magnetic-only
            // probe cell (cells:probe-non-magnetic-cells). The atom's own
            // presence in `cell` already proves the cell is non-empty, and
            // dipole::cells_field_array_*[cell] is always well-defined by
            // this point (either genuinely computed, or defaulted to 0.0 by
            // the sentinel clean-up in update_field()/hierarchical::update()
            // for a cell nobody visited), so it is always safe to copy.

            // Copy B-field from macrocell to atomistic spin
            dipole::atom_dipolar_field_array_x[atom] = dipole::cells_field_array_x[cell];
            dipole::atom_dipolar_field_array_y[atom] = dipole::cells_field_array_y[cell];
            dipole::atom_dipolar_field_array_z[atom] = dipole::cells_field_array_z[cell];

            // Unroll Hdemag field
            dipole::atom_mu0demag_field_array_x[atom] = dipole::cells_mu0Hd_field_array_x[cell];
            dipole::atom_mu0demag_field_array_y[atom] = dipole::cells_mu0Hd_field_array_y[cell];
            dipole::atom_mu0demag_field_array_z[atom] = dipole::cells_mu0Hd_field_array_z[cell];
         }

         // dipole:output-atomistic-dipole-field for the macrocell/tensor/
         // hierarchical solvers: the atomistic/atomisticfft solvers write
         // atomistic_dipole_field.txt themselves (calculate_atomistic_dipole_field(),
         // atomistic.cpp); this is the equivalent hook for everything routed
         // through this broadcast, using the one-time coordinate output added
         // to dipole::initialize() (initialize.cpp).
         if(dipole::internal::output_atomistic_dipole_field) dipole::internal::output_atomistic_dipole_fields();

         return;

      } // end of function

} // end of internal namespace

} // end of dipole namespace
