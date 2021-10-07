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

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

//------------------------------------------------------------------------
// Class function to normalise exchange interactions
//------------------------------------------------------------------------
void unitcell::exchange_template_t::normalise_exchange(){

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
         normalise_exponential_exchange();
         break;

      case internal::material_exponential:
         normalise_material_exponential_exchange();
         break;

      default:
         return;

   }

   return;

}

//------------------------------------------------------------------------
// Function to normalise exchange interactions
//------------------------------------------------------------------------
void unitcell::exchange_template_t::normalise_exponential_exchange(){

   // calculate expected sum from all nearest neighbours
   double expected_sum = 0.0;
   for(int a=0; a<ni.size(); a++){
      expected_sum += double(ni[a]);
   }

   // calculate actual sum of all interactions
   double sum = 0.0;
   for(int i=0; i<interaction.size(); i++){
      sum += interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // normalise to get same sum as for nearest neighbours
   const double inv_norm_factor = expected_sum/sum;
   for(int i=0; i<interaction.size(); i++){
      interaction[i].Jij[0][0] *= inv_norm_factor;
      interaction[i].Jij[1][1] *= inv_norm_factor;
      interaction[i].Jij[2][2] *= inv_norm_factor;
   }

   //double nsum = 0.0;
   //for(int i=0; i<interaction.size(); i++){
   //   nsum += interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   //}

   // output sum and normalised sum to screen
   //std::cout << expected_sum << "\t" << sum << "\t" << nsum << std::endl;

   return;

}

void unitcell::exchange_template_t::normalise_material_exponential_exchange(){

   int num_materials = 3;

   // Calculate expected number of interactions for different material pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double>> material_expected_sum(num_materials, std::vector<double>(num_materials, 0.0));

   for (int i = 0; i < interaction.size(); ++i){
      int min_mat = std::min(interaction[i].mat_i, interaction[i].mat_j);
      int max_mat = std::max(interaction[i].mat_i, interaction[i].mat_j);
      material_expected_sum[min_mat][max_mat] += 1.0;
   }

   // Calculate sum of interaction energies for different material category pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double>> material_sum(num_materials, std::vector<double>(num_materials, 0.0));

   for(int i=0; i<interaction.size(); ++i){
      int min_mat = std::min(interaction[i].mat_i, interaction[i].mat_j);
      int max_mat = std::max(interaction[i].mat_i, interaction[i].mat_j);
      material_sum[min_mat][max_mat] += interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // Obtain inverse normalisation factors for different material category pairs. Stored by [lowest mat cat][highest mat cat]
   std::vector <std::vector <double>> mat_inv_norm_factor(num_materials, std::vector<double>(num_materials, 0.0));
   
   for(int i = 0; i < mat_inv_norm_factor.size(); ++i){
      for (int j = 0; j < mat_inv_norm_factor.size(); ++j){
         if(j >= i){ // Only need to do this set of possibities due to symmetry of the problem
            if(material_sum[i][j] < 0.000000001) mat_inv_norm_factor[i][j] = 0; // Prevents division by zero when exchange function is set to zero
            else mat_inv_norm_factor[i][j] = material_expected_sum[i][j]/material_sum[i][j]; // Calculates the standard inv_norm_factor
         }
      }
   }
   
   // Calculate normalised exchange energies
   for(int i=0; i<interaction.size(); ++i){
      int min_mat = std::min(interaction[i].mat_i, interaction[i].mat_j);
      int max_mat = std::max(interaction[i].mat_i, interaction[i].mat_j);
      interaction[i].Jij[0][0] *= mat_inv_norm_factor[min_mat][max_mat];
      interaction[i].Jij[1][1] *= mat_inv_norm_factor[min_mat][max_mat];
      interaction[i].Jij[2][2] *= mat_inv_norm_factor[min_mat][max_mat];
   }

   return;

}



} // end if namespace unitcell
