
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
      std::vector<double> calculate_ax(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material, double unit_cell_size_x, double unit_cell_size_y, double unit_cell_size_z, std::vector <double> volume_array, std::vector <double> x_coord_array, std::vector <double> y_coord_array, std::vector <double> z_coord_array, double num_atoms_in_unit_cell)
      {
         double a2d[num_cells][num_cells] = {0.0};
         std::vector<double>  a(num_cells*num_cells,0.0);
         double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;
         for (int atom = 0; atom <num_atoms; atom++){

            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){
               if (cell_array[atom] != cell_array[neighbour_list_array[neighbour]])
               {
                  double dx = x_coord_array[atom] - x_coord_array[neighbour_list_array[neighbour]];
                  a2d[cell_array[atom]][cell_array[neighbour_list_array[neighbour]]]  =+  mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]]*dx*dx;
               //   std::cout << dx << "\t" << x_coord_array[atom] << '\t' << x_coord_array[neighbour_list_array[neighbour]]<<  "\t"<< mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]] << '\t' << a2d[cell_array[atom]][cell_array[neighbour_list_array[neighbour]]]<< std::endl;
               }
            }
         }
         int k = 0;
         for (int i =0; i < num_cells; i++)
         {
            double sum = 0;
            for (int j =0; j <num_cells; j++)
            {
            //   std::cout << a2d[i][j << '\t' << volume_array[i] << '\t' << ms[i] << '\t' << unit_cell_size_x << '\t' << atomic_volume <<std::endl;
               a[k] = (a2d[i][j]*volume_array[i])/(2*ms[i]*unit_cell_size_x*unit_cell_size_x*atomic_volume);
               k++;
            }
         //   std::cout << sum <<std::endl;
         }

         return a;
      }


      std::vector<double> calculate_ay(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material, double unit_cell_size_x, double unit_cell_size_y, double unit_cell_size_z, std::vector <double> volume_array, std::vector <double> x_coord_array, std::vector <double> y_coord_array, std::vector <double> z_coord_array, double num_atoms_in_unit_cell)
      {
         double a2d[num_cells][num_cells] = {0.0};
         std::vector<double>  a(num_cells*num_cells,0.0);
         double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;
         for (int atom = 0; atom <num_atoms; atom++){

            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){
               if (cell_array[atom] != cell_array[neighbour_list_array[neighbour]])
               {
                  double dy = y_coord_array[atom] - y_coord_array[neighbour_list_array[neighbour]];
                  a2d[cell_array[atom]][cell_array[neighbour_list_array[neighbour]]]  =+  mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]]*dy*dy;
               }
            }
         }
         int k = 0;
         for (int i =0; i < num_cells; i++)
         {
            for (int j =0; j <num_cells; j++)
            {
               a[k] = (a2d[i][j]*volume_array[i])/(2*ms[i]*unit_cell_size_y*unit_cell_size_y*atomic_volume);
               k++;
            }
         }

         return a;
      }

      std::vector<double> calculate_az(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material, double unit_cell_size_x, double unit_cell_size_y, double unit_cell_size_z, std::vector <double> volume_array, std::vector <double> x_coord_array, std::vector <double> y_coord_array, std::vector <double> z_coord_array, double num_atoms_in_unit_cell)
      {
         double a2d[num_cells][num_cells] = {0.0};
         std::vector<double>  a(num_cells*num_cells,0.0);
         double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;
         std::cout <<atomic_volume <<std::endl;
         for (int atom = 0; atom <num_atoms; atom++){

            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom]; neighbour ++){
               if (cell_array[atom] != cell_array[neighbour_list_array[neighbour]])
               {
                  double dz = z_coord_array[atom] - z_coord_array[neighbour_list_array[neighbour]];
                  a2d[cell_array[atom]][cell_array[neighbour_list_array[neighbour]]]  =+  mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]]*dz*dz;
               }
            }
         }
         int k = 0;
         for (int i =0; i < num_cells; i++)
         {
            for (int j =0; j <num_cells; j++)
            {
               a[k] = (a2d[i][j]*volume_array[i])/(2*ms[i]*unit_cell_size_z*unit_cell_size_z*atomic_volume);
               k++;
            }
         }

         return a;
      }

   }
}
