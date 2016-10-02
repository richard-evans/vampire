
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic
{

   namespace internal
   {
      std::vector<double> calculate_a(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material)
      {

         std::vector<double>  a(num_cells*num_cells,0.0);

         for (int atom = 0; atom <num_atoms; atom++){

            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){

               if (cell_array[atom] != cell_array[neighbour]) a[atom + num_materials*neighbour] = a[atom + num_materials*neighbour] +  material[type_array[atom]].Jij_matrix_SI[type_array[neighbour]];
            }
         }
         return a;
      }
   }
}
