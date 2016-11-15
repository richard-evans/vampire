
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic
{

   namespace internal
   {
      //gamma = sum(gamma)/N for each cell
      std::vector<double> calculate_gamma(int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,
                                          const std::vector<int> type_array,
                                          std::vector <mp::materials_t> material)
      {



         std::vector<double>  gamma(num_cells,0.0);
         std::vector<double>  N(num_cells,0.0);

         //gamma = 1/N sum_(i = N) gamma
         for (int atom = 0; atom <num_atoms; atom++){

            gamma[cell_array[atom]] = gamma[cell_array[atom]] + mp::material[type_array[atom]].gamma_rel;
            N[cell_array[atom]]++;
         }
         //calculates gamma/N
         for (int cell = 0; cell < num_cells; cell++){
            gamma[cell] = gamma[cell]/N[cell];
         }
         return gamma;
      }
   }
}
