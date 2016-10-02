// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"


namespace micromagnetic
{
   namespace internal
   {
      std::vector<double> calculate_tc(int num_atoms, int num_cells, std::vector<int> cell_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material)
      {

         std::vector<double>  Tc(num_cells,0.0);
         const double e = 1;
         const double kB = 1.38064852e-23;
         std::vector<double>  Jij_sum(num_cells,0.0);

         for (int atom = 0; atom <num_atoms; atom++){

            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){

               Jij_sum[cell_array[atom]] = Jij_sum[cell_array[atom]] + material[type_array[atom]].Jij_matrix_SI[type_array[neighbour]];
            }
         }



         for (int cell = 0; cell < num_cells; cell++)
         {
            Tc[cell] = Jij_sum[cell]*e/(3*kB);
         }
      return Tc;
      }
   }
}
