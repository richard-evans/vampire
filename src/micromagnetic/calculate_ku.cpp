// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic
{
   namespace internal
   {
      //ku = sum (ku) for each cell
      std::vector<double> calculate_ku(const int num_atoms,
                                       const int num_cells,
                                       std::vector<int>  cell_array,
                                       const std::vector<int> type_array,
                                       std::vector <mp::materials_t> material){

         //ku = sum ku
         std::vector<double> ku(num_cells,0.0);
         //sums over all atoms to add the ku per cell
         for (int atom = 0; atom < num_atoms; atom++)
         {
            ku[cell_array[atom]] = ku[cell_array[atom]] +  mp::material[type_array[atom]].Ku1_SI;

         }
         return ku;
      }
   }
}
