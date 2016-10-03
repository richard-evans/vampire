// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic
{

   namespace internal
   {
      std::vector<double> calculate_chi_perp(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array, std::vector <mp::materials_t> material)
      {
         std::vector<double>  chi(num_cells,0.0);

         return chi;
      }
   }
}
