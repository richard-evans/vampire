
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic{

   namespace internal{

      //calcualtes gyromagnetic ratio (gamma) for each cell from the individual atom properties.
      //gamma = sum(gamma)/N for each cell

      std::vector<double> calculate_gamma(int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                          const std::vector<int> type_array,            //1D array storing which material each atom is
                                          std::vector <mp::materials_t> material){      //class of material parameters for the atoms


         std::vector<double>  gamma(num_cells,0.0);     //1D vector storing the value of gamma for each cell
         std::vector<double>  N(num_cells,0.0);         //

         //gamma = 1/N sum_(i = N) gamma
         for (int atom = 0; atom <num_atoms; atom++){
           int cell = cell_array[atom];
           int mat = type_array[atom];
            gamma[cell] = gamma[cell] + mp::material[mat].gamma_rel;
            N[cell_array[atom]]++;
         }
         //calculates gamma/N
         for (int cell = 0; cell < num_cells; cell++){
            gamma[cell] = gamma[cell]/N[cell];
         }
         return gamma;                     //returns a 1D array of values of gamma for each cell
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
