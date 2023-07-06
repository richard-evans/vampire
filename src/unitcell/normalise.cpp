//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020, Jack B Collings 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

//------------------------------------------------------------------------
// Class function to normalise exchange interactions
//------------------------------------------------------------------------
void unitcell::exchange_template_t::normalise_exchange(std::vector < std::vector <double> > &nn_cutoff_range){

   // Select program to run
   switch(uc::internal::exchange_function){

      case internal::nearest_neighbour:
         return;
         break;

      case internal::shell:{
         return;
         break;
      }

      case internal::exponential:
         normalise_functional_exchange(nn_cutoff_range);
         break;

      case internal::material_exponential:
         normalise_functional_exchange(nn_cutoff_range);
         break;

      case internal::RKKY:
         normalise_functional_exchange(nn_cutoff_range);

      default:
         return;

   }

   return;

}

//------------------------------------------------------------------------
// Function to normalise exchange interactions
//------------------------------------------------------------------------

void unitcell::exchange_template_t::normalise_functional_exchange(std::vector < std::vector <double> > &nn_cutoff_range){

   const int num_materials = internal::material_exchange_parameters.size();

   // Calculate expected number of interactions for different material pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double> > material_expected_sum(num_materials, std::vector<double>(num_materials, 0.0));

   // Calculate sum of interaction energies for different material category pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double> > material_sum(num_materials, std::vector<double>(num_materials, 0.0));

   for (size_t i = 0; i < interaction.size(); ++i){
      int min_mat = std::min(interaction[i].mat_i, interaction[i].mat_j);
      int max_mat = std::max(interaction[i].mat_i, interaction[i].mat_j);
      if(interaction[i].rij < nn_cutoff_range[min_mat][max_mat]){
         ++material_expected_sum[min_mat][max_mat];
      }
      material_sum[min_mat][max_mat] += interaction[i].Jij[1][1]; // xx only since J is a trace anyway
   }

   // Obtain inverse normalisation factors for different material category pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double> > mat_inv_norm_factor(num_materials, std::vector<double>(num_materials, 0.0));

   for(size_t i = 0; i < mat_inv_norm_factor.size(); ++i){
      for (size_t j = 0; j < mat_inv_norm_factor.size(); ++j){
         if(j >= i){ // Only need to do this set of possibities due to symmetry of the problem
            if(fabs(material_sum[i][j]) < 0.000000001) {
               mat_inv_norm_factor[i][j] = 0; // Prevents division by zero when exchange function is set to zero
               zlog << zTs() << "unit-cell-category[" << i << "][" << j << "] interaction strength set to zero" << std::endl;
            }
            else mat_inv_norm_factor[i][j] = material_expected_sum[i][j]/material_sum[i][j]; // Calculates the standard inv_norm_factor
         }
      }
   }

   // Calculate normalised exchange energies
   for(size_t i=0; i<interaction.size(); ++i){
      int min_mat = std::min(interaction[i].mat_i, interaction[i].mat_j);
      int max_mat = std::max(interaction[i].mat_i, interaction[i].mat_j);
      interaction[i].Jij[0][0] *= mat_inv_norm_factor[min_mat][max_mat];
      interaction[i].Jij[1][1] *= mat_inv_norm_factor[min_mat][max_mat];
      interaction[i].Jij[2][2] *= mat_inv_norm_factor[min_mat][max_mat];
   }

   return;

}

} // end if namespace unitcell
