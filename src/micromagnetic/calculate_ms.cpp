//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "micromagnetic.hpp"
#include "cells.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic {
   namespace internal {

      //calculates the saturation magnetisation of each cell as a sum of the magnetic moment of each atom in the cell.
      //ms = sum muS for each cell
      std::vector<double> calculate_ms(int num_local_cells,
                                        const int num_atoms,
                                        const int num_cells,
                                        std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                        const std::vector<int> type_array,            //1D array storing which material each atom is
                                        std::vector <mp::materials_t> material,
                                        std::vector <int>local_cell_array){      //class of material parameters for the atoms

         //stores ms for each cell
         std::vector<double> ms(num_cells,0.0);
         std::vector<int> N(num_cells,0.0);
         //sums over all atoms to sum the muS per cell
         for (int atom = 0; atom < num_atoms; atom++) {
            int cell = cell_array[atom];
            int mat = type_array[atom];
            ms[cell] = ms[cell] + material[mat].mu_s_SI;
            N[cell]++;
         }
        // for (int cell = 0; cell < num_cells; cell++)
         //std::cin.get();
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &ms[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
         #endif
         return ms;           //returns a 1D array containg the saturation magnetisation of every cell
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
