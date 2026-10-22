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
#include "constants.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

// C library headers
#include <fenv.h>
#include <signal.h>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "micromagnetic.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"
#include "sim.hpp"
#include "vutil.hpp"
// dipole module headers
#include "internal.hpp"

// alias internal hierarchical namespace for brevity
namespace ha = hierarchical::internal;

namespace hierarchical{

//------------------------------------------------------------------------------
// Function to update hierarchical cell magnetization and compute dipole field
//------------------------------------------------------------------------------
void update(std::vector <double>& x_spin_array, // atomic spin directions
            std::vector <double>& y_spin_array,
            std::vector <double>& z_spin_array,
            std::vector <double>& m_spin_array, // atomic spin moment
            std::vector < bool >& magnetic){ // is magnetic

   // update hierarchical magnetization in cells
   hierarchical::internal::calculate_hierarchical_magnetisation(x_spin_array, y_spin_array, z_spin_array, m_spin_array, magnetic);

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer for the field loop
   timer.start();

   // Compute dipole fields for all cells with atoms (local cells)
   for(int lc = 0; lc < dipole::cells_num_local_cells; lc++){

      // get global cell ID from local cell list
      int cell_i = cells::cell_id_array[lc];

      // Store range of hierarchical cells contributing to dipole field for lc cell
      const int start = ha::interaction_list_start_index[lc];
      const int end = ha::interaction_list_end_index[lc];

      // Self demagnetisation factor multiplying m(i)
      const double V = dipole::cells_volume_array[cell_i];
      const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*V);
      const double self_demag = eightPI_three_cell_volume;
      //std::cout << self_demag << std::endl;

      // Cell magnetisation (Bohr magnetons) from level 0 of the hierarchical cells, whose ids are the cell ids
      const double mx_i = ha::mag_array_x[cell_i];
      const double my_i = ha::mag_array_y[cell_i];
      const double mz_i = ha::mag_array_z[cell_i];
      // std::cout << cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << std::endl;

      // Add self-demagnetisation as mu_0/4_PI * 8PI*m_cell/3V
      dipole::cells_field_array_x[cell_i] = self_demag * mx_i; //*0.0
      dipole::cells_field_array_y[cell_i] = self_demag * my_i; //*0.0
      dipole::cells_field_array_z[cell_i] = self_demag * mz_i; //*0.0

      // Loop over all cells
      for(int j = start; j<end;j++){

         // get cell ID of neighbouring cell
         int cell_j = ha::interaction_list[j];

         const double mx = ha::mag_array_x[cell_j];
         const double my = ha::mag_array_y[cell_j];
         const double mz = ha::mag_array_z[cell_j];

         //    if (cell_i == 389)  std::cout << cell_j << '\t' << mx << '\t' << my << '\t' << mz<< "\t" << cell_i << '\t' << std::endl;
         //if (cell_i == 0)std::cout<< cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << "\t" <<  cell_j << '\t' << mx << '\t' << my << '\t' << mz <<std::endl;

         // Compute dipole field contribution using dipole tensor
         dipole::cells_field_array_x[cell_i]      +=(mx*ha::rij_tensor_xx[j] + my*ha::rij_tensor_xy[j] + mz*ha::rij_tensor_xz[j]);
         dipole::cells_field_array_y[cell_i]      +=(mx*ha::rij_tensor_xy[j] + my*ha::rij_tensor_yy[j] + mz*ha::rij_tensor_yz[j]);
         dipole::cells_field_array_z[cell_i]      +=(mx*ha::rij_tensor_xz[j] + my*ha::rij_tensor_yz[j] + mz*ha::rij_tensor_zz[j]);
         // std::cout << rij_tensor_xx[j] << '\t' << rij_tensor_xy[j] << '\t' << rij_tensor_xz[j] << '\t' << rij_tensor_yy[j] << '\t' << rij_tensor_yz[j] << '\t' << rij_tensor_zz[j] << '\t' <<std::endl;

      }
      //   std::cout <<"D\t" << cell_i << '\t' << dipole::cells_field_array_x[cell_i] <<'\t' << dipole::cells_field_array_y[cell_i] <<'\t' << dipole::cells_field_array_z[cell_i] <<std::endl;

      // Multiply Hdemg by mu_0/4pi * 1e30 * mu_B to account for normalisation of magnetisation and volume in angstrom
      dipole::cells_field_array_x[cell_i] = dipole::cells_field_array_x[cell_i] * constants::dipole_prefactor;
      dipole::cells_field_array_y[cell_i] = dipole::cells_field_array_y[cell_i] * constants::dipole_prefactor;
      dipole::cells_field_array_z[cell_i] = dipole::cells_field_array_z[cell_i] * constants::dipole_prefactor;
      //std::cout << sim::time << '\t' << cell_i << '\t' << dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;
      //std::cout << "SP" << dipole::cells_field_array_x[cell_i] << "\t" << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i]<< std::endl;
      //if (cell_i == 0) std::cout << sim::time << '\t' << cell_i << '\t' <<  dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;

   }

   timer.stop();
   ha::phase_time[ha::field_loop] = timer.elapsed_time();
   timer.start();

   // Each process has the fields of the cells holding its atoms; micromagnetic cells need all cells
   if(micromagnetic::discretisation_type > 0) dipole::gather_cell_fields();

   // Copy B-field from macrocell to atomistic spins, and the demagnetising field
   // mu_0 H = B - mu_0 M (self term -1/3 mu_0 M, as in the tensor solver) for the magnetostatic energy
   for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){
      if( magnetic[atom] ){
         const int cell = cells::atom_cell_id_array[atom];
         dipole::atom_dipolar_field_array_x[atom] = dipole::cells_field_array_x[cell];
         dipole::atom_dipolar_field_array_y[atom] = dipole::cells_field_array_y[cell];
         dipole::atom_dipolar_field_array_z[atom] = dipole::cells_field_array_z[cell];
         const double mu0_m = 4.0 * M_PI * constants::dipole_prefactor / dipole::cells_volume_array[cell];
         dipole::atom_mu0demag_field_array_x[atom] = dipole::cells_field_array_x[cell] - mu0_m * ha::mag_array_x[cell];
         dipole::atom_mu0demag_field_array_y[atom] = dipole::cells_field_array_y[cell] - mu0_m * ha::mag_array_y[cell];
         dipole::atom_mu0demag_field_array_z[atom] = dipole::cells_field_array_z[cell] - mu0_m * ha::mag_array_z[cell];
      }
   }

   timer.stop();
   ha::phase_time[ha::field_distribution] = timer.elapsed_time();

   // Optionally output dipole update time to log file (commented due to call frequency of this function)
   //std::cout << "\tdone! [ " << timer.elapsed_time() << " s ]" << std::endl;
   //zlog << zTs() <<  "\tDIPOLE UPDATE. Time taken: " << timer.elapsed_time() << " s"<< std::endl;
   //std::cout << "dipole update time " << timer.elapsed_time() << " s" << std::endl;

   dipole::output_dipole_fields();

   return;

} // end of update function

//------------------------------------------------------------------------------
// Function returning the time of each phase of the last update (slowest process)
//------------------------------------------------------------------------------
std::string update_phase_times(){

   // field loop time per interaction on this process (ns), to separate load imbalance from per-interaction cost
   const double field_time = ha::phase_time[ha::field_loop];
   const double num_interactions = ha::interaction_list.size();
   const double ns_per_interaction = num_interactions > 0 ? 1.0e9 * field_time / num_interactions : 0.0;

   std::ostringstream times;
   times << " upward_s="     << vmpi::reduce_max(ha::phase_time[ha::upward_pass])
         << " reduce_s="     << vmpi::reduce_max(ha::phase_time[ha::moment_reduction])
         << " field_s="      << vmpi::reduce_max(field_time)
         << " distribute_s=" << vmpi::reduce_max(ha::phase_time[ha::field_distribution])
         << " field_min_s="  << vmpi::reduce_min(field_time)
         << " field_mean_s=" << vmpi::reduce_sum(field_time) / double(vmpi::num_processors)
         << " field_ns_per_interaction_min=" << vmpi::reduce_min(ns_per_interaction)
         << " field_ns_per_interaction_max=" << vmpi::reduce_max(ns_per_interaction);

   return times.str();

}

} // end of hierarchical namespace
