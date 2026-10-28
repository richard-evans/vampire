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

   // initialise the cells this process contributes to (the partial moments of all other cells stay zero)
   for(const int cell : ha::contributing_cells){

      ha::partial_mag_array_x[cell] = 0.0;
      ha::partial_mag_array_y[cell] = 0.0;
      ha::partial_mag_array_z[cell] = 0.0;

   }

   // If discretisation is not micromagnetic then compute cell magnetizations from atoms
   if( micromagnetic::discretisation_type != 1 ){

      //calculate total moment in each local cell looping over local atoms
      for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){

         // get cell_ID for atom
         const int cell = cells::atom_cell_id_array[atom];

         // copy spin moment to temporary variable for performance
         const double mus = m_spin_array[atom];

         // Consider only magnetic elements
         if( magnetic[atom] ){
            ha::partial_mag_array_x[cell] += x_spin_array[atom] * mus;
            ha::partial_mag_array_y[cell] += y_spin_array[atom] * mus;
            ha::partial_mag_array_z[cell] += z_spin_array[atom] * mus;
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

         ha::partial_mag_array_x[cell] = cells::mag_array_x[cell]*imuB;
         ha::partial_mag_array_y[cell] = cells::mag_array_y[cell]*imuB;
         ha::partial_mag_array_z[cell] = cells::mag_array_z[cell]*imuB;

      }

   }

   //--------------------------------------------------------------------------------------
   // loop over the hierarchical levels in order, computing partial cell magnetizations at
   // level L for the cells holding local cells (partial sums in other cells are zero)
   //--------------------------------------------------------------------------------------
   for (const int cell : ha::upward_pass_cells){

      // determine which lower level cells L-1 are in each higher level cell
      int start_cell_in_cell = ha::cells_in_cells_start_index[cell];
      int end_cell_in_cell = ha::cells_in_cells_end_index[cell];

      // loop over all L-1 level cells and accumulate magnetization
      for (int cell_in_cell = start_cell_in_cell; cell_in_cell < end_cell_in_cell; cell_in_cell++){

         // get local cell ID of subcell
         int subcell = ha::cells_in_cells[cell_in_cell];

         // accumulate cell magnetizations up the hierarchical tree
         ha::partial_mag_array_x[cell] += ha::partial_mag_array_x[subcell];
         ha::partial_mag_array_y[cell] += ha::partial_mag_array_y[subcell];
         ha::partial_mag_array_z[cell] += ha::partial_mag_array_z[subcell];
      //   std::cout << ha::partial_mag_array_z[cell] <<std::endl;
      }
      //  std::cout << ha::partial_mag_array_x[cell] << '\t' << ha::partial_mag_array_y[cell] << '\t' << ha::partial_mag_array_z[cell] << '\t' << std::endl;

   } // end of cell loop

   timer.stop();
   ha::phase_time[ha::upward_pass] = timer.elapsed_time();
   timer.start();

   // Finally sum the contributions of all processors into the shared cell moments
   const size_t num_contributing_cells = ha::contributing_cells.size();
   for(size_t i = 0; i < num_contributing_cells; i++){
      const int cell = ha::contributing_cells[i];
      ha::mag_reduction_values[i]                              = ha::partial_mag_array_x[cell];
      ha::mag_reduction_values[i +     num_contributing_cells] = ha::partial_mag_array_y[cell];
      ha::mag_reduction_values[i + 2 * num_contributing_cells] = ha::partial_mag_array_z[cell];
   }
   ha::mag_reduction.sum(ha::mag_reduction_values, {&ha::mag_array_x, &ha::mag_array_y, &ha::mag_array_z});

   timer.stop();
   ha::phase_time[ha::moment_reduction] = timer.elapsed_time();

   return;

}

} // end of internal namespace
} // end of hierarchical namespace
