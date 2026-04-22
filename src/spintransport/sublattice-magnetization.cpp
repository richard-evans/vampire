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

   //---------------------------------------------------------------------------
   // reset magnetization vector to zero
   //---------------------------------------------------------------------------
   for( int sl = 0; sl < num_sublattices; sl++ ){
      std::fill(st::internal::cell_sl_magnetization[sl].begin(), st::internal::cell_sl_magnetization[sl].end(), 0.0);
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
      st::internal::cell_sl_magnetization[sl][3*cell+0] += mm*atoms_x_spin_array[atom];
      st::internal::cell_sl_magnetization[sl][3*cell+1] += mm*atoms_y_spin_array[atom];
      st::internal::cell_sl_magnetization[sl][3*cell+2] += mm*atoms_z_spin_array[atom];

   }

   // for(int cell = 0; cell < st::internal::cell_sl_magnetization[0].size()/3; cell++){
   //    for( int sl = 0; sl < num_sublattices; sl++ ){
   //       std::cout << "cell: " << cell << "\t" <<
   //                 st::internal::cell_sl_magnetization[sl][3*cell+0] << "\t" <<
   //                 st::internal::cell_sl_magnetization[sl][3*cell+1] << "\t" <<
   //                 st::internal::cell_sl_magnetization[sl][3*cell+2] << std::endl;
   //    }
   // }

   //---------------------------------------------------------------------------
   // Reduce cell material magnetizations on all processors
   //---------------------------------------------------------------------------
   #ifdef MPICF
      // cast to int for MPI
      int bufsize = 3*st::internal::total_num_cells;
      for( int sl = 0; sl < num_sublattices; sl++ ){
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_magnetization[sl][0], bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
   #endif

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
