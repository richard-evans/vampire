
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic
{

   namespace internal
   {
      std::vector<double> calculate_gamma(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array, std::vector <mp::materials_t> material)
      {

         std::vector<double>  gamma(num_cells*num_cells,0.0);

         for (int atom = 0; atom <num_atoms; atom++){

            gamma[cell_array[atom]] = gamma[cell_array[atom]] + material[type_array[atom]].temperature_rescaling_alpha;
         }
         return gamma;
      }
   }
}
