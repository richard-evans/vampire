//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "spintransport.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{
namespace internal{

//---------------------------------------------------------------------------------------------------------
// Function to update cell material magnetizations m / m_s^0
//---------------------------------------------------------------------------------------------------------
void calculate_cell_sublattice_magnetization(const unsigned int num_local_atoms,             // number of local atoms
                                             const std::vector<double>& atoms_x_spin_array,  // x-spin vector of atoms
                                             const std::vector<double>& atoms_y_spin_array,  // y-spin vector of atoms
                                             const std::vector<double>& atoms_z_spin_array,  // z-spin-vector of atoms
                                             const std::vector<double>& atoms_m_spin_array   // moment of atoms
){

   // define local constants for number of cells sublattices to avoid repeated access to global variables
   const int num_sublattices = st::internal::num_sublattices;
   const int num_cells = st::internal::total_num_cells;

   //---------------------------------------------------------------------------
   // reset magnetization vector to zero
   //---------------------------------------------------------------------------
   for( int cell = 0; cell < num_cells; cell++ ){
      for( int sl = 0; sl < num_sublattices; sl++ ){
         st::internal::cell_sl_magnetization_x[cell][sl] = 0.0;
         st::internal::cell_sl_magnetization_y[cell][sl] = 0.0;
         st::internal::cell_sl_magnetization_z[cell][sl] = 0.0;
      }
   }

   //---------------------------------------------------------------------------
   // loop over all atoms and determine cell magnetizations (can OpenMP)
   //---------------------------------------------------------------------------
   for(unsigned int atom = 0; atom < num_local_atoms; atom++){

      // get cell id
      const uint64_t cell = st::internal::atom_in_cell[atom];

      // get magnetic moment (muB)
      const double mm = atoms_m_spin_array[atom];

      // get sublattice of atom
      const int sl = st::internal::atom_sublattice[atom];

      // add magnetization to cell
      st::internal::cell_sl_magnetization_x[cell][sl] += mm*atoms_x_spin_array[atom];
      st::internal::cell_sl_magnetization_y[cell][sl] += mm*atoms_y_spin_array[atom];
      st::internal::cell_sl_magnetization_z[cell][sl] += mm*atoms_z_spin_array[atom];

   }

   // for(int cell = 0; cell < num_cells; cell++){
   //    for( int sl = 0; sl < num_sublattices; sl++ ){
   //       std::cout << "cell: " << cell << "\t" <<
   //                 st::internal::cell_sl_magnetization_x[cell][sl] << "\t" <<
   //                 st::internal::cell_sl_magnetization_y[cell][sl] << "\t" <<
   //                 st::internal::cell_sl_magnetization_z[cell][sl] << std::endl;
   //    }
   // }

   //---------------------------------------------------------------------------
   // Reduce cell material magnetizations on all processors
   //---------------------------------------------------------------------------
   #ifdef MPICF
      // cast to int for MPI
      int bufsize = num_sublattices;
      for( int cell = 0; cell < num_cells; cell++ ){
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_magnetization_x[cell][0], bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_magnetization_y[cell][0], bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_magnetization_z[cell][0], bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
   #endif

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
