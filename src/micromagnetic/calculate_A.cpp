
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"
#include "create.hpp"
#include <stdlib.h>
namespace micromagnetic
{

   namespace internal
   {

      //---------------------------------------------------------------------------------
      //
      //                       A = 1/4*sum((Jij) *(x_i - x_j)^2)
      //This calculates the interaction between neighbouring cells therefore only interactions where i and j
      //are in different macrocells are summed.
      //
      // we take some components of the exchange field in here and calulate them as part of a therefore we return
      //
      //      A = (1/2)*sum((Jij)(x_i - x_j)^2) * (V_cell/V_atomic) * (1/l_cell) * (1/Ms)
      //
      //-----------------------------------------------------------------------------

      std::vector< double > calculate_a(int num_atoms,
                                        int num_cells,
                                        int num_materials,
                                        std::vector<int> cell_array,
                                        std::vector<int> neighbour_list_array,
                                        std::vector<int> neighbour_list_start_index,
                                        std::vector<int> neighbour_list_end_index,
                                        const std::vector<int> type_array,
                                        std::vector <mp::materials_t> material,
                                        double unit_cell_size,
                                        std::vector <double> volume_array,
                                        std::vector <double> x_coord_array,
                                        std::vector <double> y_coord_array,
                                        std::vector <double> z_coord_array,
                                        double num_atoms_in_unit_cell){

         std::vector< std::vector< double> > a2d(num_cells, num_cells, 0.0); 
         std::vector<double> a;
         a.resize(num_cells*num_cells,0.0);

         //calculates the atomic volume  = volume of one cell/number of atoms in a unitcell = atomic volume
         const double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;

         //loops over all atoms and all neighbours and checks whether i and j are in the same cell or neighbouring cells
         for (int atom = 0; atom <num_atoms; atom++){
            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom] +1; neighbour ++){
               if (cell_array[atom] != cell_array[neighbour_list_array[neighbour]]){
                  //if atom i and j are not in the same cell calculate A
                  double dx = x_coord_array[atom] - x_coord_array[neighbour_list_array[neighbour]];
                  double dy = y_coord_array[atom] - y_coord_array[neighbour_list_array[neighbour]];
                  double dz = z_coord_array[atom] - z_coord_array[neighbour_list_array[neighbour]];
                  double d = sqrt(dx*dx + dy*dy + dz*dz);
                  a2d[cell_array[atom]][cell_array[neighbour_list_array[neighbour]]] += mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]]*d*d;
                }
            }
         }
         k = 0;
         // loops over all cells to turn the 2D array into a 1D array
         // multiplys A by cell size/2Ms*V_Atomic to ad din the terms of H_Ex
         for (int i =0; i < num_cells; i++){
            double cell_size = pow(volume_array[i],1./3.);
            for (int j =0; j <num_cells; j++){
               if ((ms[j] != 0) && (ms[i] != 0) && (a2d[i][j] < 0.5)){
                  k = j*num_cells +i;
                  a[k] = (a2d[i][j]*cell_size)/(2*ms[i]*atomic_volume);
               }
               k++;
            }
         }

         int k1,k2;
         k1 = 0;
         //sums over all interactions to check interaction between cell i j = interaction cell j i
         //non symetric interactions not realistic
         for (int i = 0; i < num_cells; i ++){
            for (int j = 0; j < num_cells; j++){
               k2 = (j*num_cells) + i;
               k1 = (i*num_cells) + j;
               if (a2d[i][j] != a2d[j][i]) std::cout << "error: Non symetric exchange" <<"\t"  <<  i << '\t' << j << "\t"  <<  a2d[i][j]<<"\t"  <<  a2d[j][i] <<std::endl;
            }
         }
         return a;
      }
   }
}
