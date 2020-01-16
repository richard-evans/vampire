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

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic{

   namespace internal{

      // Calculates the intrinsic damping constant of each cell from the atomsitic parameters in each cell.

      //alpha = sum(alpha)/N for each cell

      std::vector<double> calculate_alpha(int num_local_cells,
                                          int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,                 //1D array storing which cell each atom is in
                                          const std::vector<int> type_array,           //1D array storing which material each atom is
                                          std::vector <mp::materials_t> material,
                                          std::vector<int> local_cell_array){     //class of material parameters for the atoms

         std::vector<double>  alpha(num_cells,0.0);                                    //vectors to store alpha
         std::vector<double>  N(num_cells,0.0);                                        //vector stores the number of atoms per micromagnetic cell

         //sums over all atoms to calulcate the sum of alpha and number of atoms for each cell
         for (int atom = 0; atom <num_atoms; atom++){
            int cell = cell_array[atom];
            int mat = type_array[atom];
            alpha[cell] += mp::material[mat].alpha;
            N[cell] += 1.0;
         }

         // Reduce sum of alpha and N on all processors
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &alpha[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &N[0],     num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         // calculates the average alpha per cell
         for (int cell = 0; cell < num_cells; cell++){
            alpha[cell] = alpha[cell]/N[cell];
         }

         // reduce final alpha values on all cells
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &alpha[0], num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         #endif

         return alpha;          //return an array of damping constants for each cell

      }

   } //closes the internal namspace
}  //closes the micromagnetic namespace
