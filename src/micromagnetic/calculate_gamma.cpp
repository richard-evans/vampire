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

      //calculates gyromagnetic ratio (gamma) for each cell from the individual atom properties.
      //gamma = sum(gamma)/N for each cell

      std::vector<double> calculate_gamma(int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                          const std::vector<int> type_array,            //1D array storing which material each atom is
                                          std::vector <mp::materials_t> material,
                                          int num_local_cells,
                                          std::vector<int>local_cell_array){      //class of material parameters for the atoms


         std::vector<double>  gamma(num_cells,0.0);     //1D vector storing the value of gamma for each cell
         std::vector<double>  N(num_cells,0.0);         //

         //gamma = 1/N sum_(i = N) gamma
         for (int atom = 0; atom <num_atoms; atom++){
            int cell = cell_array[atom];
            int mat = type_array[atom];
            gamma[cell] = gamma[cell] + mp::material[mat].gamma_rel;
            N[cell_array[atom]]++;
         }

         // Reduce sum of gamma and N on all processors
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &gamma[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &N[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         //calculates gamma/N
         for (int cell = 0; cell < num_cells; cell++){
            gamma[cell] = gamma[cell]/N[cell];
         }

         // reduce final gamma values on all cells
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &gamma[0], num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         #endif

         return gamma;                     //returns a 1D array of values of gamma for each cell

      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
