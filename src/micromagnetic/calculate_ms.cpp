
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic
{
   namespace internal
   {
      //ms = sum muS for each cell
      std::vector<double> calculate_ms(const int num_atoms, const int num_cells, std::vector<int> cell_array, const std::vector<int> type_array,  std::vector <mp::materials_t> material)
      {
         //stores ms for each cell
         std::vector<double> ms(num_cells,0.0);
         //sums over all atoms to sum the muS per cell
         for (int atom = 0; atom < num_atoms; atom++)
         {
            ms[cell_array[atom]] = ms[cell_array[atom]] + material[type_array[atom]].mu_s_SI;
         }
      //   std::cout << "\t" << ms[0] <<std::endl;
         return ms;
      }
   }
}
