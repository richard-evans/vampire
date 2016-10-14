
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic
{

   namespace internal
   {
      std::vector<double> calculate_tc(int num_atoms, int num_cells,std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material)
      {

         std::vector<double>  J(num_cells,0.0);
         std::vector<double>  N(num_cells,0.0);
         std::vector<double>  Tc(num_cells,0.0);
         const double e = 0.86;
         const double kB = 1.38064852e-23;

         for (int atom = 0; atom <num_atoms; atom++){
            N[cell_array[atom]]++;
            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){
               J[cell_array[atom]] = J[cell_array[atom]] +  mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]];
            }
         }

         for (int cell = 0; cell < num_cells; cell++)
         {
            Tc[cell] = -J[cell]*e/(3*kB*N[cell]);

         }
         return Tc;
      }
   }
}
