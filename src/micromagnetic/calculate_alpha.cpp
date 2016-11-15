
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic
{

   namespace internal
   {
      //alpha = sum(alpha)/N for each cell
      std::vector<double> calculate_alpha(int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,
                                          const std::vector<int> type_array,
                                          std::vector <mp::materials_t> material){
         //stores alpha per cell
         std::vector<double>  alpha(num_cells,0.0);
         //stores the number of atoms in each cell
         std::vector<double>  N(num_cells,0.0);

         //sums over all atoms to calulcate the sum of alpha for each cell
         for (int atom = 0; atom <num_atoms; atom++){
            alpha[cell_array[atom]] = alpha[cell_array[atom]] + mp::material[type_array[atom]].alpha;
            N[cell_array[atom]]++;
         }

         //divides the sum of alpha by N
         for (int cell = 0; cell < num_cells; cell++){
            alpha[cell] = alpha[cell]/N[cell];
         }
         return alpha;
      }
   }
}
