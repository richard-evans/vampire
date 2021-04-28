

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "hierarchical.hpp"

// exchange module headers
#include "internal.hpp"

namespace hierarchical{

   //------------------------------------------------------------------------------
   // Functions to return unrolled dipole tensor
   //------------------------------------------------------------------------------
   std::vector<double> get_tensor_1D_xx(){
      return hierarchical::internal::rij_tensor_xx;
   }

   std::vector<double> get_tensor_1D_xy(){
      return hierarchical::internal::rij_tensor_xy;
   }

   std::vector<double> get_tensor_1D_xz(){
      return hierarchical::internal::rij_tensor_xz;
   }
   std::vector<double> get_tensor_1D_yy(){
      return hierarchical::internal::rij_tensor_yy;
   }
   std::vector<double> get_tensor_1D_yz(){
      return hierarchical::internal::rij_tensor_yz;
   }
   std::vector<double> get_tensor_1D_zz(){
      return hierarchical::internal::rij_tensor_zz;
   }   
    std::vector<int> get_interaction_list(){
      return hierarchical::internal::interaction_list;
   }   
    std::vector<int> get_interaction_list_start(){
      return hierarchical::internal::interaction_list_start_index;
   } 
    std::vector<int> get_interaction_list_end(){
      return hierarchical::internal::interaction_list_end_index;
   } 
    std::vector<double> get_cell_mag_x(){
      return hierarchical::internal::mag_array_x;
   } 
    std::vector<double> get_cell_mag_y(){
      return hierarchical::internal::mag_array_y;
   } 
    std::vector<double> get_cell_mag_z(){
      return hierarchical::internal::mag_array_z;
   } 
    int get_num_cells(){
      return hierarchical::internal::total_num_cells;
   } 

} // end of dipole namespace




