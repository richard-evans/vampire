//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins, Andrea Meo and Richard F L Evans 2019.
//   All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include "constants.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>

// Vampire headers
#include "cells.hpp" // needed for dp::cell_id_array but to be removed
#include "vmpi.hpp"
#include "vutil.hpp"
#include "atoms.hpp"

// hierarchical module headers
#include "hierarchical.hpp"
#include "internal.hpp"
#include "micromagnetic.hpp"

// alias internal hierarchical namespace for brevity
namespace ha = hierarchical::internal;

namespace hierarchical{
namespace internal{

//-----------------------------------------------------------------------------
// Function for calculate magnetisation in hierarchical cells
//-----------------------------------------------------------------------------
// As cells are duplicated on all processors, but local cells are not unique,
// the accumulation of moments must be done directly from the atomic positions.
// Each process accumulates its contribution (partial moments), and the
// contributions are summed into the cell moments shared by each region.
//-----------------------------------------------------------------------------
void calculate_hierarchical_magnetisation(std::vector <double>& x_spin_array, // atomic spin directions
                                          std::vector <double>& y_spin_array,
                                          std::vector <double>& z_spin_array,
                                          std::vector <double>& m_spin_array, // atomic spin moment
                                          std::vector < bool >& magnetic){ // is magnetic

   // time the local accumulation and upward pass, then the reduction
   vutil::vtimer_t timer;
   timer.start();

   // Partial moments are compact: index i holds this process's contribution to cell ha::contributing_cells[i]
   // (initialise.cpp). Initialise them; all other cells receive nothing from this process.
   const size_t num_contributing_cells = ha::contributing_cells.size();
   std::fill(ha::partial_mag_array_x.begin(), ha::partial_mag_array_x.end(), 0.0);
   std::fill(ha::partial_mag_array_y.begin(), ha::partial_mag_array_y.end(), 0.0);
   std::fill(ha::partial_mag_array_z.begin(), ha::partial_mag_array_z.end(), 0.0);

   // If discretisation is not micromagnetic then compute cell magnetizations from atoms
   if( micromagnetic::discretisation_type != 1 ){

      //calculate total moment in each local cell looping over local atoms
      for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){

         // Consider only magnetic elements
         if( magnetic[atom] ){

            // compact index of the atom's cell
            const int i = ha::atom_contributing_index[atom];

            // copy spin moment to temporary variable for performance
            const double mus = m_spin_array[atom];

            ha::partial_mag_array_x[i] += x_spin_array[atom] * mus;
            ha::partial_mag_array_y[i] += y_spin_array[atom] * mus;
            ha::partial_mag_array_z[i] += z_spin_array[atom] * mus;
         }

      } // End of atom loop

   }
   // Otherwise for micromagnetic simulations use cell arrays
   else {

      // inverse Bohr magneton
      const double imuB = constants::i_muB;

      // initialise locally integrated cells to cell magnetization values
      for (int lc = 0; lc < micromagnetic::number_of_micromagnetic_cells; lc++){

         const int cell = micromagnetic::list_of_micromagnetic_cells[lc];
         const int i = ha::micromagnetic_contributing_index[lc];

         ha::partial_mag_array_x[i] = cells::mag_array_x[cell]*imuB;
         ha::partial_mag_array_y[i] = cells::mag_array_y[cell]*imuB;
         ha::partial_mag_array_z[i] = cells::mag_array_z[cell]*imuB;

      }

   }

   //--------------------------------------------------------------------------------------
   // loop over the hierarchical levels in order, computing partial cell magnetizations at
   // level L for the cells holding local cells, from their children that hold local cells
   // (the partial sums of the other children are zero)
   //--------------------------------------------------------------------------------------
   const size_t num_upward_pass_cells = ha::upward_pass_index.size();
   for (size_t u = 0; u < num_upward_pass_cells; u++){

      const int i = ha::upward_pass_index[u];

      // accumulate cell magnetizations up the hierarchical tree
      for (int c = ha::upward_children_start[u]; c < ha::upward_children_start[u+1]; c++){
         const int child = ha::upward_children[c];
         ha::partial_mag_array_x[i] += ha::partial_mag_array_x[child];
         ha::partial_mag_array_y[i] += ha::partial_mag_array_y[child];
         ha::partial_mag_array_z[i] += ha::partial_mag_array_z[child];
      }

   } // end of cell loop

   timer.stop();
   ha::phase_time[ha::upward_pass] = timer.elapsed_time();
   timer.start();

   // Finally sum the contributions of all processors into the shared cell moments
   for(size_t i = 0; i < num_contributing_cells; i++){
      ha::mag_reduction_values[i]                              = ha::partial_mag_array_x[i];
      ha::mag_reduction_values[i +     num_contributing_cells] = ha::partial_mag_array_y[i];
      ha::mag_reduction_values[i + 2 * num_contributing_cells] = ha::partial_mag_array_z[i];
   }
   if(ha::moment_exchange == ha::moment_exchange_t::point_to_point) ha::exchange_moments(); // exchange.cpp
   else ha::mag_reduction.sum(ha::mag_reduction_values, {&ha::mag_array_x, &ha::mag_array_y, &ha::mag_array_z});

   timer.stop();
   ha::phase_time[ha::moment_reduction] = timer.elapsed_time();

   return;

}

} // end of internal namespace
} // end of hierarchical namespace
