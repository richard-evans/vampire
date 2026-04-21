//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spintransport.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{
namespace internal{

//---------------------------------------------------------------------------------------------------------
// Function to update cell magnetizations m / m_s^0
//---------------------------------------------------------------------------------------------------------
void calculate_cell_magnetization(const unsigned int num_local_atoms,            // number of local atoms
                                  const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
                                  const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
                                  const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
                                  const std::vector<double>& atoms_m_spin_array  // moment of atoms
                               ){


   //---------------------------------------------------------------------------
   // reset magnetization vector to zero
   //---------------------------------------------------------------------------
   std::fill(st::internal::cell_magnetization.begin(), st::internal::cell_magnetization.end(), 0.0);

   //---------------------------------------------------------------------------
   // loop over all atoms and determine cell magnetizations (can OpenMP)
   //---------------------------------------------------------------------------
   for(unsigned int atom = 0; atom < num_local_atoms; atom++){

      // get cell id
      const uint64_t cell = st::internal::atom_in_cell[atom];

      // get magnetic moment (muB)
      const double mm = atoms_m_spin_array[atom];

      // add magnetization to cell
      st::internal::cell_magnetization[3*cell+0] += mm*atoms_x_spin_array[atom];
      st::internal::cell_magnetization[3*cell+1] += mm*atoms_y_spin_array[atom];
      st::internal::cell_magnetization[3*cell+2] += mm*atoms_z_spin_array[atom];

   }

   //---------------------------------------------------------------------------
   // Reduce cell magnetizations on all processors
   //---------------------------------------------------------------------------
   #ifdef MPICF
      // cast to int for MPI
      int bufsize = 3*st::internal::total_num_cells;
      MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_magnetization[0], bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
