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
#include <algorithm>
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
   const int num_sublattices = spin_transport::internal::num_sublattices;

   //---------------------------------------------------------------------------
   // reset magnetization vector to zero
   //---------------------------------------------------------------------------
   std::fill(spin_transport::internal::cell_sl_magnetization_x.begin(), spin_transport::internal::cell_sl_magnetization_x.end(), 0.0);
   std::fill(spin_transport::internal::cell_sl_magnetization_y.begin(), spin_transport::internal::cell_sl_magnetization_y.end(), 0.0);
   std::fill(spin_transport::internal::cell_sl_magnetization_z.begin(), spin_transport::internal::cell_sl_magnetization_z.end(), 0.0);

   //---------------------------------------------------------------------------
   // loop over all atoms and determine cell magnetizations (can OpenMP)
   //---------------------------------------------------------------------------
   for(unsigned int atom = 0; atom < num_local_atoms; atom++){

      // get cell id
      const uint64_t cell = spin_transport::internal::atom_in_cell[atom];

      // get magnetic moment (muB)
      const double mm = atoms_m_spin_array[atom];

      // get sublattice of atom
      const int sl = spin_transport::internal::atom_sublattice[atom];

      // flat index: sl is the fastest-varying index
      const int idx = (int)cell * num_sublattices + sl;

      // add magnetization to cell
      spin_transport::internal::cell_sl_magnetization_x[idx] += mm*atoms_x_spin_array[atom];
      spin_transport::internal::cell_sl_magnetization_y[idx] += mm*atoms_y_spin_array[atom];
      spin_transport::internal::cell_sl_magnetization_z[idx] += mm*atoms_z_spin_array[atom];

   }

   // for(int cell = 0; cell < num_cells; cell++){
   //    for( int sl = 0; sl < num_sublattices; sl++ ){
   //       std::cout << "cell: " << cell << "\t" <<
   //                 spin_transport::internal::cell_sl_magnetization_x[cell][sl] << "\t" <<
   //                 spin_transport::internal::cell_sl_magnetization_y[cell][sl] << "\t" <<
   //                 spin_transport::internal::cell_sl_magnetization_z[cell][sl] << std::endl;
   //    }
   // }

   //---------------------------------------------------------------------------
   // Reduce cell material magnetizations on all processors
   //---------------------------------------------------------------------------
   #ifdef MPICF
      const int num_cells = spin_transport::internal::total_num_cells;
      const int num_cells_x_sl = num_cells * num_sublattices;
      // single reduction per component over the entire flat array
      MPI_Allreduce(MPI_IN_PLACE, &spin_transport::internal::cell_sl_magnetization_x[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &spin_transport::internal::cell_sl_magnetization_y[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &spin_transport::internal::cell_sl_magnetization_z[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
